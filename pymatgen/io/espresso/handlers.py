# coding: utf-8

from __future__ import unicode_literals, division

"""
"""


__author__ = "David Cossey"
__version__ = "0.1"
__maintainer__ = "David Cossey"
__email__ = "dcossey014@gmail.com"
__status__ = "Beta"
__date__ = "10/17/17"

import math
import logging
import os, glob

from fireworks import LaunchPad, Firework
from custodian.custodian import ErrorHandler
from custodian.utils import backup
from pymatgen.io.espresso.inputs import PWInput

logger = logging.getLogger(__name__)

QE_BACKUP_FILES  = ['in', 'out', 'pp_in', 'pp_out']

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
                #print("line: {}".format(line))
                logger.debug("line: {}".format(line))
                l = line.strip()
                for err, msgs in QuantumEspressoErrorHandler.error_msgs.items():
                    for msg in msgs:
                        if l.find(msg) != -1:
                            #print("\n\nfound error: {}\nat line: {}".format(msg, l))
                            logger.debug("\n\nfound error: {}\nat line: {}".format(msg, l))
                            self.errors.add(err)
            #print("Found these errors: {}".format(self.errors))
            logger.info("Found these errors: {}".format(self.errors))
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

