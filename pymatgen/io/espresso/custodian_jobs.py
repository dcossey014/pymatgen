from __future__ import division, unicode_literals

__author__ = 'David Cossey, Gary Kedziora'
__version__ = '0.1'
__maintainer__ = 'Gary Kedziora'
__email__ = 'dcossey014@gmail.com; gary.kedziora@engilitycorp.com'
__date__ = '10/17/17'


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
from pymatgen.io.espresso.inputs import PWInput, PwInputError
from pymatgen.io.espresso.outputs import EspressoRun, PWOutput
from pymatgen.io.bgw.outputs import BgwRun, BgwParserError
from fireworks import Firework, FireTaskBase, FWAction, explicit_serialize, Workflow, LaunchPad
from custodian.custodian import Job, Custodian

def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)

logger = logging.getLogger(__name__)


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
        os.environ['ESPRESSO_NPROCS'] = self.pw_cmd[2]
        pw_exe = [i for i in self.pw_cmd if 'pw.x' in i][0]
        iotk_cmd = os.path.join(os.path.dirname(pw_exe), 'iotk')
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

