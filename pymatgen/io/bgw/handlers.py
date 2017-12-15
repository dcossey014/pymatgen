# coding: utf-8

from __future__ import unicode_literals, division

"""
"""


__author__ = "David Cossey"
__version__ = "0.2"
__maintainer__ = "David Cossey"
__email__ = "dcossey014@gmail.com"
__status__ = "Beta"
__date__ = "7/01/17"

import math, datetime
import logging
import os, glob

from fireworks import LaunchPad, Firework
from custodian.custodian import ErrorHandler
from custodian.utils import backup
from pymatgen.io.bgw.inputs import BgwInput
from pymatgen.io.bgw.interfaces import BgwFirework, BgwCustodianTask
from pymatgen.io.bgw.interpreter import BgwModder
from custodian.ansible.interpreter import Modder

logger = logging.getLogger(__name__)

BGW_BACKUP_FILES = ["OUT.*", "*.inp", "FW*", "bgw.log", "pbs.log"]

class BgwErrorHandler(ErrorHandler):
    """
    Error handler for BerkeleyGW Jobs. 
    """

    is_monitor = True

    error_msgs = {
                'error':    ['MPT ERROR', 'ERROR'
                            ]
                }

    #TODO: Fix runtype with glob
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

        input_files = glob.glob("*.inp")
        if len(input_files) > 1:
            logger.info("found more than one input file.  Using {}".format(
                input_files[0]))
            print("found more than one input file.  Using {}".format(
                input_files[0]))
        elif len(input_files) < 1:
            logger.info("Could not find an input file.")
            print("Could not find an input file.")
            raise NameError("Could not find an input file.")

        self.run_type = input_files[0].split('.')[0]
        self.bgwi = BgwInput.from_file(input_files[0])
        self.output_filename = "OUT.{}".format(self.run_type[:3])

    def check(self):
        '''
        Checks output file for errors.
        '''

        self.errors = set()
        with open(self.output_filename, 'r') as fout:
            for line in fout:
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
        logger = logging.getLogger(__name__)

        input_files = glob.glob("*.inp")
        if len(input_files) > 1:
            logger.info("found more than one input file.  Using {}".format(
                input_files[0]))
            print("found more than one input file.  Using {}".format(
                input_files[0]))
        elif len(input_files) < 1:
            logger.info("Could not find an input file.")
            print("Could not find an input file.")
            raise NameError("Could not find an input file.")

        self.run_type = input_files[0].split('.')[0]
        logger.debug("init for BgwMemoryHandler\nrun_type: {}".format(self.run_type))
        self.bgwi = BgwInput.from_file(input_files[0])
        self.errors = []
        self.output_filename = "OUT.{}".format(self.run_type[:3])

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


class WalltimeErrorHandler(ErrorHandler):
    """
    Check if a run is nearing the walltime. If so, terminate the running process
    so that the BgwCustodianTask() can create a new Firework to continue the
    job if necessary.  You can specify the walltime either in the init (
    which is unfortunately necessary for SGE and SLURM systems. If you happen
    to be running on a PBS system and the PBS_WALLTIME variable is in the run
    environment, the wall time will be automatically determined if not set. If 
    using AFRL PBS templates, this Environmental Variable should be generated
    automatically from the PBS script.
    """
    is_monitor = True
    is_terminating = True

    # This handler will be unrecoverable, but custodian shouldn't raise an
    # error
    raises_runtime_error = False

    def __init__(self, walltime=None, buffer_time=300):
        """
        Initializes the handler with a buffer time.

        Args:
            walltime (int): Total walltime in seconds. If this is None and
                the job is running on a PBS system, the handler will attempt to
                determine the walltime from the PBS_WALLTIME environment
                variable. If the wall time cannot be determined or is not
                set, this handler will have no effect.
            buffer_time (int): The min amount of buffer time in secs at the
                end that the STOPCAR will be written. The STOPCAR is written
                when the time remaining is < the higher of 3 x the average
                time for each ionic step and the buffer time. Defaults to
                300 secs, which is the default polling time of Custodian.
                This is typically sufficient for the current ionic step to
                complete. But if other operations are being performed after
                the run has stopped, the buffer time may need to be increased
                accordingly.
        """
        if walltime is not None:
            self.walltime = walltime

        elif "PBS_WALLTIME" in os.environ:
            self.walltime = os.environ["PBS_WALLTIME"]
            wt = self.walltime.split(":")
            wt_delta = datetime.timedelta(hours=int(wt[0]), 
                    minutes=int(wt[1]), 
                    seconds=int(wt[2])) if len(wt) > 1 else None
            self.walltime = wt_delta.total_seconds() if wt_delta else int(self.walltime)

        else:
            self.walltime = None

        self.buffer_time = buffer_time
        self.start_time = datetime.datetime.now()
        self.prev_check_time = self.start_time
        ##DEC DEBUG
        print("printing: Found this walltime: {}".format(self.walltime))
        logger.info("logging: Found this walltime: {}".format(self.walltime))
        ##END DEBUG

    def check(self):
        if self.walltime:
            run_time = datetime.datetime.now() - self.start_time
            total_secs = run_time.total_seconds()

            # If the remaining time is less than average time for 3 ionic
            # steps or buffer_time.
            time_left = self.walltime - total_secs
            if time_left < self.buffer_time:
                return True
        
        ##DEC DEBUG
        print("Time left = {}".format(time_left))
        logger.info("Time left = {}".format(time_left))
        ##END DEBUG

        return False

    def correct(self):
        return {"errors": ["Walltime reached"], "actions": None}

