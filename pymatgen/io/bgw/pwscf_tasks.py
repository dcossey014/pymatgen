from __future__ import division, unicode_literals

__author__ = 'David Cossey, Gary Kedziora'
__version__ = '0.0'
__maintainer__ = 'Gary Kedziora'
__email__ = 'dcossey014@gmail.com; gary.kedziora@engilitycorp.com'
__date__ = '3/09/16'


import os
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
from pymongo import MongoClient

from pymatgen import Structure
from pymatgen.io.pwscf import PWInput, PWInputError, PWOutput
from pymatgen.io.bgw.outputs import BgwRun, EspressoRun, BgwParserError
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
        print "in BGWJob.__init__"
        logger.debug("in BGWJob.__init__")

        self.bgw_cmd=bgw_cmd
        self.pp_cmd = pp_cmd
        self.output_file=output_file
        print("bgw_cmd: {}".format(self.bgw_cmd))
        logger.debug("bgw_cmd: {}".format(self.bgw_cmd))

    def setup(self):
        print "in BGWJob.setup"

    def run(self):
        cmd=list(self.bgw_cmd)
        print "in BGWJob.run, cmd = ", cmd
        logger.debug("in BGWJob.run, cmd = ", cmd)
        with open(self.output_file, 'w') as f:
            p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        return p

    def postprocess(self):
        print "in BGWJob.postprocess"
        logger.debug("in BGWJob.postprocess")
        if self.pp_cmd:
            #print("pp_cmd: {}".format(self.pp_cmd))
            #print("type of pp_cmd: {}".format(type(self.pp_cmd)))
            logger.debug("pp_cmd: {}".format(self.pp_cmd))
            logger.debug("type of pp_cmd: {}".format(type(self.pp_cmd)))
            cmd = self.pp_cmd.split() if isinstance(self.pp_cmd,
                        (str, unicode)) else list(self.pp_cmd)
            #print "running PostProcessing with: ", cmd
            logger.debug("running PostProcessing with: ", cmd)
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
    required_params = ['bgw_cmd', 'fout']
    optional_params = ['pp_cmd', 'mpi_cmd', 'handlers', 'handler_params']

    def run_task(self, fw_spec):
        FORMAT = "%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - \n\t\t%(message)s\n"
        logging.basicConfig(filename="bgw.log", level=logging.DEBUG, format=FORMAT)
        logger = logging.getLogger(__name__)

        mpix = self.get('mpi_cmd', 'mpiexec_mpt -n 36').split()
        bgwx = list(mpix)
        bgwx.extend([self.get('bgw_cmd')])
        fout = self.get('fout')
        ppx = self.get('pp_cmd', None)
        prev_dirs = fw_spec.get('PREV_DIRS', {})

        #print('prev dirs: {}'.format(prev_dirs))

        #Set up Handlers and their Parameters
        handler_param_dict = {'BgwMemoryHandler': ['run_type']}

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
    
            print("runnung {} Handler with params: {}".format(n, np))
            logger.info("runnung {} Handler with params: {}".format(n, np))

        logger.info("creating BGWJob")
        job = BGWJob(bgwx, pp_cmd=ppx, output_file=fout)
        c = Custodian(handlers=handlers, validators=[], jobs=[job], monitor_freq=3)
        output = c.run()
        print("output returned from Custodian: {}".format(output))
        return FWAction(stored_data=output[0], mod_spec=[{'_set': {
                #'PREV_DIRS': prev_dirs,
                'PREV_DIRS->BGW->{}'.format(os.path.basename(self.get(
                                'bgw_cmd')).split('.')[0]): os.getcwd()}}])


@explicit_serialize
class BgwDB(FireTaskBase):
    required_params = ['config_file']
    optional_params = ['insert_to_db']

    def run_task(self, fw_spec):
        self.db_config = loadfn(os.path.join(os.environ['HOME'], self.get('config_file')))
        self.upload = self.get('insert_to_db', False)
        self.database = self.db_config.get('database', 'BGW_DATA')
        self.collection = self.db_config.get('collection', 'DEFAULT')
        self.username = self.db_config.get('username', None)
        self.password = self.db_config.get('password', None)
        self.ssl_ca_file = self.db_config.get('ssl_ca_file', None)


        self.prev_dirs = fw_spec.get("PREV_DIRS", None)
        self.esp_dir = os.path.dirname(self.prev_dirs.get("ESPRESSO", {}).get("scf", None))
        self.bgw_dirs = self.prev_dirs.get("BGW", {})

        if self.esp_dir:
            self.esp_data = EspressoRun(self.esp_dir)
        if self.bgw_dirs:
            for i in self.bgw_dirs:
                setattr(self, i, BgwRun(self.bgw_dirs[i]))

        # Upload to MongoDB or return a PrettyPrint Dictionary
        if self.upload:
            self.insert_db(self.as_dict())
        else:
            #pp = pprint.PrettyPrinter(indent=2)
            #pp.pprint(self.as_dict())
            pass


    def get_dict(self, esp_dir, bgw_dirs):
        esp_data = EspressoRun(esp_dir)
        d = {}
        for i,r in enumerate(bgw_dirs):
            out_file = glob.glob(os.path.join(bgw_dirs[r], "OUT.*"))

            if len(out_file) > 1:
                raise BgwParserError(
                        "Found more than one output file for "
                        "Runtype: {}".format(tmp_out.runtype),
                        {'err': 'Duplicate Run', 'file': out_file[-1]} )

            tmp_out = BgwRun(out_file[0])
            
            if tmp_out.runtype not in d.keys() and tmp_out.timings:
                d.update(tmp_out.as_dict())

            elif tmp_out.runtype in d.keys():
                raise BgwParserError(
                        "Found more than one output file for "
                        "Runtype: {}".format(tmp_out.runtype),
                        {'err': 'Duplicate Run', 'file': r} )
            else:
                raise BgwParserError(
                        "Could not parse timings.  Make sure the run "
                        "completed successfully.",
                        {'err': 'Incomplete Run'} )
        
        return (d, esp_data.as_dict())


    def insert_db(self, run_data):
        connection = MongoClient(self.db_config['host'], self.db_config['port'],
                        ssl_ca_certs=self.ssl_ca_file)
        
        db = connection[self.database]
        if self.username:
            db_auth = connection['admin']
            db_auth.authenticate(self.username, self.password)

        collection = db[self.collection]
        collection.insert_one(run_data)

    def as_dict(self):
        esp_run = self.esp_data.as_dict()
        d = {
            'Structure': esp_run['structure'],
            'ESPRESSO': esp_run,
            'BGW': {}
            }

        for i in self.bgw_dirs:
            bgw_run = getattr(self, i)
            d['BGW'][i] = bgw_run.as_dict()
        return d


@explicit_serialize
class WritePwscfInputTask(FireTaskBase):

    required_params = ["structure", "pseudo", "control"]
    optional_params = ["system", "electrons", "ions", "cell", "kpoints_mode", 
                        "kpoints_grid", "kpoints_shift"]

    def run_task(self, fw_spec):

        #prev_dir = fw_spec.get('PREV_DIR', None)

        if isinstance(self["structure"], Structure):
            s = self["structure"]
        elif isinstance(self["structure"], dict):
            s = Structure.from_dict(self["structure"])
        else:
            s = Structure.from_file(os.path.join(prev_dir, self["structure"]))

        p=self.get("pseudo")
        c=self.get("control")
        sy=self.get("system", None)
        e=self.get("electrons", None)
        ions=self.get("ions", None)
        cell=self.get("cell", None)
        kmode=self.get("kpoints_mode", "automatic")
        kgrid=self.get("kpoints_grid", (1,1,1))
        kshift=self.get("kpoints_shift", (0,0,0))

        pwin=PWInput(s,control=c,pseudo=p,system=sy,electrons=e,ions=ions,cell=cell,kpoints_mode=kmode,
                        kpoints_grid=kgrid,kpoints_shift=kshift)

        pwin.write_file("in")


@explicit_serialize
class SimplePWTask(FireTaskBase):

    def run_task(self, fw_spec):
        print "in SimplePWTask.run_task"
        # gk: this is a site dependent script.  It would be better to use 
        # the Quantum Espress executable directly.
        job=PWJob(["pw.x"])
        c=Custodian(handlers=[], validators=[], jobs=[job])
        output=c.run()


class PWJob(Job):
    def __init__(self, pw_cmd, pw2bgw_cmd=None, output_file="out"):
        print "in PWJob.__init__"
        self.pw_cmd=pw_cmd
        self.pw2bgw_cmd = pw2bgw_cmd
        self.output_file=output_file
        print("pw_cmd: {}".format(self.pw_cmd))

    def setup(self):
        print "in PWJob.setup"

    def run(self):
        cmd=list(self.pw_cmd)
        print "in PWJob.run, cmd = ", cmd
        fin=open("in",'r')
        with open(self.output_file, 'w') as f:
            p = subprocess.Popen(cmd, stdin=fin, stdout=f, stderr=subprocess.STDOUT)
        return p

    def postprocess(self):
        print "in PWJob.postprocess"
        print "converting charge-density.dat into XML format"
        iotk_cmd = os.path.join(os.path.dirname(self.pw_cmd[-1]), 'iotk')
        print("iotk cmd: {}".format(iotk_cmd))
        chg_dat = glob.glob('*.save/charge-density.dat')[0]
        chg_xml = os.path.join(os.path.dirname(chg_dat), 'charge-density.xml')
        p = subprocess.Popen([iotk_cmd, 'convert', chg_dat, chg_xml])
        p.wait()

        if self.pw2bgw_cmd:
            cmd=list(self.pw2bgw_cmd)
            print "running PostProcessing with: ", cmd
            fin=open('pp_in', 'r')
            with open('pp_out', 'w') as f:
                p = subprocess.Popen(cmd, stdin=fin, stdout=f, stderr=subprocess.STDOUT)
            p.wait()

