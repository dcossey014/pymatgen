# coding: utf-8

from __future__ import unicode_literals, division

"""
"""


__author__ = "David Cossey"
__version__ = "0.1"
__maintainer__ = "David Cossey"
__email__ = "dcossey014@gmail.com"
__status__ = "Beta"
__date__ = "7/01/17"

import math
import logging
import os

from fireworks import LaunchPad, Firework
from custodian.custodian import ErrorHandler
from custodian.utils import backup
from pymatgen.io.pwscf import PWInput
from pymatgen.io.bgw.inputs import BgwInput
from pymatgen.io.bgw.interfaces import BgwFirework, BgwCustodianTask
from pymatgen.io.bgw.interpreter import BgwModder
from custodian.ansible.interpreter import Modder

logger = logging.getLogger(__name__)

QE_BACKUP_FILES  = ['in', 'out', 'pp_in', 'pp_out']
BGW_BACKUP_FILES = ["OUT.*", "*.inp", "FW*", "bgw.log", "pbs.log"]

class QuantumEspressoErrorHandler(ErrorHandler):
    '''
    Error handler for Quantum Espresso Jobs
    '''

    is_monitor=True
    error_msgs = {
            'ndiag':[
                'Error in routine check_para_diag',
                'Too few bands for required ndiag'
                ],
            'n_plane_waves':[
                'Error in routine n_plane_waves',
                'No plane waves found: running on too many processors?'
                ],
            'g-vectors': [
                'Error in routine  cdiaghg (155)'
                ],
            'cholesky': [
                'Error in routine  cdiaghg (727)',
                'problems computing cholesky'
                ],
            'error': [
                'MPT ERROR', 'ERROR'
                ]
            }

    def __init__(self):
        """
        Initializes with an output file name.

        Args:
            output_filename (str): This is the file where the stdout for nwchem
                is being redirected. The error messages that are checked are
                present in the stdout. Defaults to "mol.nwout", which is the
                default redirect used by :class:`custodian.nwchem.jobs
                .NwchemJob`.
        """
        self.pwi = PWInput.from_file("in")
        self.output_filename = "out"

    def check(self):
        '''
        Checks output file for errors.
        '''

        self.errors = set()
        with open(self.output_filename, 'r') as fout:
            print("opened '{}' file for reading".format(self.output_filename))
            import os
            print("PWD: {}".format(os.path.abspath('.')))
            for line in fout:
                print("line: {}".format(line))
                l = line.strip()
                for err, msgs in QuantumEspressoErrorHandler.error_msgs.items():
                    for msg in msgs:
                        if l.find(msg) != -1:
                            print("\n\nfound error: {}\nat line: {}".format(msg, l))
                            self.errors.add(err)
            print("Found these errors: {}".format(self.errors))
        return len(self.errors) > 0

    def correct(self):
        backup(QE_BACKUP_FILES)
        actions = []

        ''' # Example
        if self.errors.intersection(['ndiag']):
            actions.append({'dict': 'INCAR', 
                        "action": {'_set': {'key': 'val'}}})
        '''
        return {"errors": list(self.errors), "actions": actions}

    def __str__(self):
        return "QuantumEspressoErrorHandler"

class BgwErrorHandler(ErrorHandler):
    """
    Error handler for BerkeleyGW Jobs. 
    """

    is_monitor = True

    error_msgs = {
                'error':    ['MPT ERROR', 'ERROR'
                            ]
                }

    def __init__(self, run_type=""):
        """
        Initializes with an output file name.

        Args:
            output_filename (str): This is the file where the stdout for nwchem
                is being redirected. The error messages that are checked are
                present in the stdout. Defaults to "mol.nwout", which is the
                default redirect used by :class:`custodian.nwchem.jobs
                .NwchemJob`.
        """
        self.run_type = runtype
        self.bgwi = BGWInput.from_file("{}.inp".format(run_type))
        self.output_filename = "OUT.{}".format(run_type[:3])

    def check(self):
        '''
        Checks output file for errors.
        '''

        self.errors = set()
        with open(self.output_filename, 'r') as fout:
            for line in f:
                l = line.strip()
                for err, msgs in BgwErrorHandler.error_msgs.items():
                    for msg in msgs:
                        if l.find(msg) != -1:
                            self.errors.add(err)
        return len(self.errors) > 0

    def correct(self):
        backup(BGW_BACKUP_FILES)
        actions = []

        return {"errors": self.errors, "actions": actions}

    def __str__(self):
        return "BgwErrorHandler"


class BgwMemoryHandler(ErrorHandler):
    """
    Error handler that verifies BerkeleyGW Jobs have the Memory Resources necessary to complete. 
    """

    is_monitor = True
    is_terminating = True

    def __init__(self, run_type="absorption"):
        """
        Initializes with an output file name.

        Args:
            output_filename (str): This is the file where the stdout for nwchem
                is being redirected. The error messages that are checked are
                present in the stdout. Defaults to "mol.nwout", which is the
                default redirect used by :class:`custodian.nwchem.jobs
                .NwchemJob`.
        """
        logger = logging.getLogger(__name__)
        #logger.debug("init for BgwMemoryHandler\nrun_type: {}".format(run_type))

        self.run_type = run_type
        self.bgwi = BgwInput.from_file("{}.inp".format(run_type))
        self.errors = []
        self.output_filename = "OUT.{}".format(run_type[:3])

    def check(self):
        '''
        Checks output file for errors.
        '''

        # Reset Memory Required value
        self.mem_req = 0.0

        #logger.debug("checking for memory in handler")
        with open(self.output_filename, 'r') as fout:
            for line in fout:
                l = line.strip()
                if l.find("MB per PE") != -1:
                    if l.find("Memory available:") != -1:
                        self.mem_avail = float(l.split()[-4])
                        logger.debug('setting available memory: mem_avail: {}'.format(
                                        self.mem_avail))

                    else:
                        logger.debug('')
                        self.mem_req += float(l.split()[-4])
                        logger.debug('setting memory required: {}'.format(self.mem_req))

                elif l.find('Running with') != -1 and l.find('MPI task') != -1:
                    self.ncpus = int(l.split()[2])
                    logger.debug('setting # cpus: \t ncpus: {}'.format(self.ncpus))
                elif l.find("Started") != -1:
                    break
        
        logger.info("\n\navail mem: {}\nrequired mem: {}\n\n".format(
                    self.mem_avail, self.mem_req))

        if self.mem_req > self.mem_avail:
            self.errors.append("Insufficient Memory")
        
        return len(self.errors) > 0


    def correct(self):
        backup(BGW_BACKUP_FILES)
        action = []
        actions = []
        self.lp = LaunchPad.from_file(os.path.join(os.environ["HOME"],
                                        ".fireworks", "my_launchpad.yaml"))
        self.Firework = Firework.from_file('./FW.json')
        self.fw_id = self.Firework.fw_id
        self.queue_spec = self.Firework.spec.get('_queueadapter', {})
        logger.info("FW_ID: {}\tQueue_spec: {}".format(self.fw_id, self.queue_spec))
        for n,i in enumerate(self.Firework.tasks):
            if isinstance(i, BgwCustodianTask): 
                mpi_cmd = self.Firework.tasks[n]['mpi_cmd']

        ppnode = int(self.queue_spec.get('ppnode',
                        int(os.environ.get('BC_MPI_TASKS_ALLOC', self.ncpus)) / 
                        int(os.environ.get('BC_NODE_ALLOC', 1))))

        mem_node_avail = int(os.environ.get('BC_MEM_PER_NODE', 
                                        self.mem_avail * ppnode))

        logger.info("ppnode: {}\nmem_node_avail: {}".format(ppnode, mem_node_avail))

        if mem_node_avail == 0:
            raise Exception("Unable to determine Memory per node available.")

        self.nnodes = int(math.ceil((self.mem_req*self.ncpus)/mem_node_avail))
        nprocs = int(self.nnodes * ppnode)
        new_mpi = mpi_cmd.split()[:-1]
        new_mpi.append(str(nprocs))
        new_mpi = ' '.join(new_mpi)

        logger.debug("Number of nodes needed: {}".format(self.nnodes))
        logger.debug("setting new mpi_cmd to: {}".format(new_mpi))
        self.queue_spec.update({'nnodes': self.nnodes})
        logger.debug("queue_spec after update: {}".format(self.queue_spec))

        # Set the new Queue Adapter parameters
        actions.append({'launchpad': self.lp, 'fw_id': self.fw_id,
                'action': {"_queueadapter": self.queue_spec}}) 
        action.append({"_set": {"_queueadapter": self.queue_spec}})

        # Set to launch in the same directory as failed run
        actions.append({'launchpad': self.lp, 'fw_id': self.fw_id,
                'action': {"_launch_dir": os.getcwd()}})
        action.append({"_set": {"_launch_dir": os.getcwd()}})

        # Change mpi_cmd used for calculation to reflect new value
        for i,t in enumerate(self.Firework.tasks):
            if isinstance(t, BgwCustodianTask):
                actions.append({"launchpad": self.lp, 'fw_id': self.fw_id,
                    'action': {"_tasks.{}.mpi_cmd".format(i): new_mpi}})
                action.append({"_set": {"_tasks.{}.mpi_cmd".format(i): new_mpi}})
                    
        BgwModder().apply_actions(actions)
        return {"errors": self.errors, "actions": action}

    def __str__(self):
        return "BgwMemoryErrorHandler"
