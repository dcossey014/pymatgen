from __future__ import division, unicode_literals

__author__ = 'David Cossey, Gary Kedziora'
__version__ = '0.2'
__maintainer__ = 'Gary Kedziora'
__email__ = 'dcossey014@gmail.com; gary.kedziora@engilitycorp.com'
__date__ = '10/17/17'


import os, re
import sys 
import six 
import glob
import shutil
import string
import logging
import datetime
import traceback
import pprint

import subprocess
import shlex

from monty.serialization import loadfn
from monty.json import json
from pymongo import MongoClient

from pymatgen import Structure
from pymatgen.io.bgw.outputs import BgwRun, BgwParserError
from pymatgen.io.espresso.outputs import EspressoRun
from fireworks import Firework, FireTaskBase, FWAction, explicit_serialize, Workflow, LaunchPad
from custodian.custodian import Job, Custodian

def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)

logger = logging.getLogger(__name__)

class BGWJob(Job):
    '''
    Docstring
    '''
    def __init__(self, bgw_cmd, pp_cmd=None, output_file="out"):
        logger.debug("in BGWJob.__init__")

        self.bgw_cmd=bgw_cmd
        self.pp_cmd = pp_cmd
        self.output_file=output_file
        logger.debug("bgw_cmd: {}".format(self.bgw_cmd))

    def setup(self):
        print "in BGWJob.setup"

    def run(self):
        cmd=list(self.bgw_cmd)
        logger.debug("in BGWJob.run, cmd = {}".format(cmd))
        with open(self.output_file, 'w') as f:
            p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        return p

    def postprocess(self):
        logger.debug("in BGWJob.postprocess")
        if self.pp_cmd:
            #print("pp_cmd: {}".format(self.pp_cmd))
            #print("type of pp_cmd: {}".format(type(self.pp_cmd)))
            logger.debug("pp_cmd: {}".format(self.pp_cmd))
            logger.debug("type of pp_cmd: {}".format(type(self.pp_cmd)))
            cmd = self.pp_cmd.split() if isinstance(self.pp_cmd,
                        (str, unicode)) else list(self.pp_cmd)
            #print "running PostProcessing with: ", cmd
            logger.debug("running PostProcessing with: {}".format(cmd))
            with open("pp.out", 'w') as f:                
                p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
            p.wait()


@explicit_serialize
class BgwCustodianTask(FireTaskBase):
    '''
    Custodian Task for Running BerkeleyGW calculations.

    Required Parameters:
        bgw_cmd (str):  Executable for the BerkeleyGW program to be executed.
        fout (str):     Output filename to store results

    Optional Parameters:
        pp_cmd (str):   Executable for any BerkeleyGW post processing commands.
        mpi_cmd (str):  String of the MPI command used to execute BGW.  i.e.
                        "mpiexec_mpt -n 36"
        handlers (list) List of handlers to be used with the BGW Job
        hander_params (dict)    Dictionary of the parameters to be used with
                                the handlers list above
    '''
    required_params = ['bgw_cmd', 'mpi_cmd', 'fout'] 
    optional_params = ['pp_cmd', 'handlers', 'handler_params']

    def run_task(self, fw_spec):
        FORMAT = "%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - \n\t\t%(message)s\n"
        logging.basicConfig(filename="bgw.log", level=logging.DEBUG, format=FORMAT)
        logger = logging.getLogger(__name__)

        mpix = self.get('mpi_cmd').split()
        bgwx = list(mpix)
        bgwx.extend([self.get('bgw_cmd')])
        fout = self.get('fout')
        ppx = self.get('pp_cmd', None)
        prev_dirs = fw_spec.get('PREV_DIRS', {})

        #print('prev dirs: {}'.format(prev_dirs))

        #Set up Handlers and their Parameters
        handler_param_dict = {'BgwMemoryHandler': ['run_type'], 
                'WalltimeErrorHandler': ['walltime', 'buffer_time'],
                'BgwErrorHandler': ['run_type']}

        hnames = self.get('handlers', [])
        handler_params = self.get('handler_params', {})

        print("handler names: {}".format(hnames))
        print("handler_params: {}".format(handler_params))

        logger.info("handler names: {}".format(hnames))
        logger.info("handler params: {}".format(handler_params))
        handlers = []

        for n in hnames:
            logger.debug("working with {} handler".format(n))
            np = {}
            for m in handler_params:
                if m in handler_param_dict[n]:
                    np[m] = handler_params[m]

            logger.debug("found these parameters for {} handler:  {}".format(
                            n, np))

            if np: 
                handlers.append(load_class("pymatgen.io.bgw.handlers", n)(**np))
            else:
                handlers.append(load_class("pymatgen.io.bgw.handlers", n)())
    
            print("running {} Handler with params: {}".format(n, np))
            logger.info("runnung {} Handler with params: {}".format(n, np))

        logger.info("creating BGWJob")
        job = BGWJob(bgwx, pp_cmd=ppx, output_file=fout)
        c = Custodian(handlers=handlers, validators=[], jobs=[job], monitor_freq=3)
        output = c.run()
        print("output returned from Custodian: {}".format(output))

        #TODO: Check this code to make sure it works.
        # Attempt to prevent multiple epsilon/sigma jobs from overwriting each
        # other
        fw = Firework.from_file('FW.json')
        name = fw.name
        m = re.match("[\w\s]+(\d)$", name)
        if m:
            n_val = "-{}".format(m.group(1))
        else:
            n_val = ""

        return FWAction(stored_data=output[0], mod_spec=[{'_set': {
                #'PREV_DIRS': prev_dirs,
                'PREV_DIRS->BGW->{}{}'.format(os.path.basename(self.get(
                                'bgw_cmd')).split('.')[0], n_val): 
                                os.getcwd()}}])


@explicit_serialize
class BgwDB(FireTaskBase):
    required_params = ['config_file']
    optional_params = ['insert_to_db']

    def run_task(self, fw_spec):
        self.db_config = loadfn(os.path.join(os.environ['HOME'], self.get('config_file')))
        self.upload = self.get('insert_to_db', False)
        self.database = self.db_config.get('database', 'BGW_DATA')
        self.collection = self.db_config.get('collection', 'calculations')
        self.username = self.db_config.get('username', None)
        self.password = self.db_config.get('password', None)
        self.ssl_ca_file = self.db_config.get('ssl_ca_file', None)


        self.prev_dirs = fw_spec.get("PREV_DIRS", None)
        self.esp_dir = os.path.dirname(os.path.dirname(self.prev_dirs.get("ESPRESSO", 
                        {}).get("scf", None)))
        self.bgw_dirs = self.prev_dirs.get("BGW", {})

        # Call parsers and put them into class attributes (dictionaries)

        # This is where the Espresso output processing gets called
        if self.esp_dir:
            self.esp_data = EspressoRun(self.esp_dir)

        # This is where the BGW output processing gets called
        if self.bgw_dirs:
            for i in self.bgw_dirs:
                setattr(self, i, BgwRun(i,self.bgw_dirs[i]))

        # Upload to MongoDB or return a PrettyPrint Dictionary
        if self.upload:
            self.insert_db(self.get_dict())
        else:
            pp = pprint.PrettyPrinter(indent=2)
            pp.pprint(self.get_dict())


    def insert_db(self, run_data):
        connection = MongoClient(self.db_config['host'], self.db_config['port'],
                        ssl_ca_certs=self.ssl_ca_file)
        
        db = connection[self.database]
        if self.username:
            db_auth = connection['admin']
            db_auth.authenticate(self.username, self.password)


        # Separate dictionaries into separate collections.
        bgw_dict = run_data.pop('BGW')
        esp_dict = run_data.pop('ESPRESSO')
        abs_tabs = {'eigenvalues': bgw_dict['absorption']['Output'].pop('eigenvalues')}
        abs_tabs['eigenvalues no eh'] = bgw_dict['absorption']['Output'].pop('eigenvalues no eh')
        abs_tabs['absorbtion bandstructure on fine grid'] = bgw_dict['absorption']['Output'].pop(
                                                'absorbtion bandstructure on fine grid')

        # Insert Top Level data into 'calculations' Database
        collection = db[self.collection]
        ob_id = collection.insert_one(run_data)

        # Insert reference id of Top Level Data to BGW and ESPRESSO 
        # Dictionaries for References and Queries
        for i in [bgw_dict, esp_dict, abs_tabs]:
            i[u'ref_id'] = ob_id.inserted_id

        # Upload each dictionary to their respective collection
        collection = db['berkeleygw']
        print("inserting bgw_dict")
        #print("bgw_dict: {}".format(bgw_dict))
        collection.insert_one(bgw_dict)
        
        collection = db['espresso']
        print("inserting esp_dict")
        #print("esp_dict: {}".format(esp_dict))
        collection.insert_one(esp_dict)

        collection = db['absorption_tables']
        print("inserting absorption tables")
        #print("absorption tables: {}".format(abs_tabs))
        collection.insert_one(abs_tabs)


    def get_dict(self):
        esp_run = self.esp_data.as_dict()
        with open(os.path.join(self.esp_dir, 'FW.json'), 'r') as fin:
            data = json.load(fin)

        kpts = {}
        for task in esp_run['input']['mf_tasks']:
            kpts[task] = len(esp_run['band_data']['kpoints'][task])

        d = {
            'Structure': esp_run['structure'],
            'ESPRESSO': esp_run,
            'pretty_formula' : self.esp_data.pretty_formula,
            'spacegroup'  : self.esp_data.spacegroup[0],
            'spacegroup_number' : self.esp_data.spacegroup[1],
            'kpoints_fine' : self.esp_data.kpts_fine,
            'kpoints_coarse' : self.esp_data.kpts_coarse,
            'number of kpoints': kpts,
            'BGW': {},
            'started_on'   : data['created_on'],
            'completed_on' : datetime.datetime.utcnow().isoformat(),
            'directories' : {'ESPRESSO': self.esp_dir, 'BGW': self.bgw_dirs}
            }

        for i in self.bgw_dirs:
            bgw_run = getattr(self, i)
            d['BGW'][i] = bgw_run.as_dict()

            d['BGW'][i]['number of kpoints'] = {}
            if 'epsilon' in i:
                d['BGW'][i]['number of kpoints']['wfn'] = kpts['wfn']
                d['BGW'][i]['number of kpoints']['wfnq'] = kpts['wfnq']

            if 'sigma' in i or 'kernel' in i or 'absorption' in i:
                d['BGW'][i]['number of kpoints']['wfn_co'] = kpts['wfn_co']

            if 'absorption' in i:
                d['BGW'][i]['number of kpoints']['wfn_fi'] = kpts['wfn_fi']
                d['BGW'][i]['number of kpoints']['wfnq_fi'] = kpts['wfnq_fi']

        return d


