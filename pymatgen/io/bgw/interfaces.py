import os, sys
import glob
import shutil
from math import ceil

import subprocess
from pymatgen.io.bgw.kgrid import generate_kpath
from monty.serialization import loadfn
from fireworks import Firework, FireTaskBase, FWAction, \
                        explicit_serialize, FileTransferTask, Workflow, LaunchPad
from custodian import Custodian
from pymatgen import Structure
from pymatgen.io.bgw.inputs import BgwInput, BgwInputTask
from pymatgen.io.bgw.custodian_jobs import BGWJob, BgwDB, BgwCustodianTask

def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)

class BgwFirework():
    '''
    This is an interface to Fireworks for use with BGW Workflows.

    Required Params:
        -bgw_task (obj):    FireTask Object used to create a Firework.
                            2 FireTasks are automatically created upon
                            constructing the BgwFirework object from a
                            BgwInput object: [<BgwInputTask>, 
                            <BgwCustodianTask>].

    Optional Params:
        -bgw_cmd (str):     Path to BGW Executable to use for the calculation
        -mpi_cmd (str):     MPI command to use for the BGW Custodian Task.
                            Default: 'mpiexec_mpt -n 36'
        -name (str):        Name given to FireWork. Default: Bgw FW

gk:      ppx is not needed, since sigma writes eqp0.dat and eqp1.dat automatically

        -ppx (str):         Post-Processing command.  This is needed for 
                            Sigma calculations to post process for absorption
                            calculations.  Default: None.
        -handlers (list):   List of Handlers to use for BGW Custodian Task.
                            Default: None
        -handler_params (dict): Dictionary of Handler Parameters to use with 
                                the Handlers. Default: {}
        -config_file (str): Patch to configuration file. Default: None
        -complex (bool):    Use complex executables for calculations.
                            Default: False
                            
    '''
    
    def __init__(self, bgw_task, name="Bgw FW", bgw_cmd=None,
                handlers=['BgwErrorHandler', 'BgwMemoryHandler',
                'WalltimeErrorHandler'], mpi_cmd=None,
                handler_params=None, ppx=None,
                config_file=None, complex=False):
        #TODO: Fix BGW_cmd and handlers.  Correct OUT file configuration.
        if config_file:
            config_dict = loadfn(config_file)
        elif os.path.exists(os.path.join(os.environ['HOME'], 
                                    'bgw_interface_defaults.yaml')):
            config_dict = loadfn(os.path.join(os.environ['HOME'],
                                    'bgw_interface_defaults.yaml'))
        else:
            config_dict = {}

        self.name = name
        self.handlers=handlers if handlers else []
        self.handler_params=handler_params if handler_params else {}
        runtype = bgw_task.params['out_file'].split('/')[-1].split('.')[0]
        self.out_file = "OUT.{}".format(runtype[:3])
        self.bgw_cmd = bgw_cmd
        self.ppx = ppx 
        if mpi_cmd:
            self.mpi_cmd = mpi_cmd
        else:
            print("Error: No mpi_cmd given. "
                    "Please input an mpi_cmd for the job.  \n"
                    "Usage: bgw_fw = BgwFirework(<bgw_task>, mpi_cmd: <mpi_cmd>, "
                    "name=<name>, complex='complex|real', config_file='config_file', "
                    "handlers=[handlers], handler_params={key: value}")
            sys.exit(1)

        # check paramers before job submission to not waste time on queue
        bgw_task.check_params()

        if config_file:
            config_dict = loadfn(config_file)
        elif os.path.exists(os.path.join(os.environ['HOME'], 
                                    'bgw_interface_defaults.yaml')):
            config_dict = loadfn(os.path.join(os.environ['HOME'],
                                    'bgw_interface_defaults.yaml'))
        else:
            config_dict = {}

        if config_dict:
            self.custodian_opts = config_dict.get('CUSTODIAN_PARAMS', {})
            if self.custodian_opts.get('handlers', []):
                self.handlers.extend(self.custodian_opts.get('handlers', []))
            self.handler_params.update(self.custodian_opts.get('handler_params', {}))

            bgw_dir = config_dict.get('BGW_DIR', '')
            if not self.bgw_cmd and bgw_dir:
                if complex:
                    self.bgw_cmd = os.path.join(bgw_dir, "{}.{}".format( 
                        bgw_task.params['run_type'], 'cplx.x'))
                else:
                    self.bgw_cmd = os.path.join(bgw_dir, "{}.{}".format(
                        bgw_task.params['run_type'], 'real.x'))

        self.tasks=[bgw_task,BgwCustodianTask(handlers=self.handlers,           
                handler_params=self.handler_params, bgw_cmd=self.bgw_cmd,
                fout=self.out_file, mpi_cmd=self.mpi_cmd, 
                pp_cmd=self.ppx)] if isinstance(bgw_task, 
                                BgwInputTask) else [bgw_task]
        self.Firework = Firework(self.tasks, name=self.name)

        # Try to establish connection with LaunchPad
        try:
            self.LaunchPad = LaunchPad.from_file(os.path.join(os.environ["HOME"], 
                                            ".fireworks", "my_launchpad.yaml"))
        except:
            self.LaunchPad = None

    def add_task(self, task):
        '''
        Function used to add another FireTask to the Firework.  If
        given task is a BgwInput Object, it will create the BgwInputTask
        as well as add another BgwCustodianTask() to the FireTask list.
        '''

        if isinstance(task, BgwInputTask):
            self.tasks.extend([task, BgwCustodianTask(handlers=self.handlers,
                    handler_params=self.handler_params)])
        else:
            self.tasks.append(task)
        self.Fireworks = Firework(self.tasks, name=self.name)


    def to_file(self, filename):
        self.Firework.to_file(filename)


    def add_fw_to_launchpad(self):
        if self.LaunchPad:
            self.Workflow=Workflow([self.Firework])
            self.LaunchPad.add_wf(self.Workflow)
        else:
            print("No connection to LaunchPad. \n"
                    "Use 'to_file(<filename>)' to write a yaml file\n"
                    "to manually add Firework to LaunchPad later.\n")

    def add_spec(self, key, val):
        # To use nested dictionary entries, use "->"
        # as the separator for nested dictionaries
        # i.e. '_queueadapter->walltime' as the key 
        # to set {'_queueadapter': {'walltime': 'val'}}
        def get_nested_dict(input_dict, key):
            current = input_dict
            toks = key.split("->")
            n = len(toks)
            for i, tok in enumerate(toks):
                if tok not in current and i < n - 1:
                    current[tok] = {}
                elif i == n - 1:
                    return current, toks[-1]
                current = current[tok]

        def get_nested_dict2(input_dict, val):
            current = input_dict
            for k, v in val.items():
                if isinstance(v, dict):
                    if k not in current.keys():
                        current[k] = {}
                    current = current[k]
                    get_nested_dict2(current, v)
                else:
                    current[k] = v            

        if isinstance(val, dict):
            if key not in self.Firework.spec:
                self.Firework.spec[key] = {}
            current = self.Firework.spec[key]
            get_nested_dict2(current, val)            

        else:
            (d, k) = get_nested_dict(self.Firework.spec, key)
            d[k] = val

    def add_handler(self, handler, **kwargs):
        '''
        Member function to add handler and handler options to 
        all Fireworks in the BgwFirework.tasks list.

        Example (assuming fw1 is the firework object):
        fw1.add_handler('WalltimeHandler', wall_time=3600, buffer_time=300)
        fw1.add_handler('FrozenJobErrorHandler')
        '''
        for i,j in enumerate(self.Firework.tasks):
            if isinstance(j, BgwCustodianTask):
                cust_handlers = j.get('handlers', [])
                cust_params = j.get('handler_params', {})
                if handler not in cust_handlers:
                    j['handlers'].append(handler)
                    j['handler_params'].update(kwargs)
        self.Firework.spec['_tasks']=[t.to_dict() for t in
                self.Firework.tasks] #re-initialize FW.spec


class BgwWorkflow():
    '''
    A BgwWorkflow encapsulates multiple BgwFirework objects into a Single Workflow.
    If the kwarg "dependency" is not set, it will create a Sequential Workflow where 
    the next Firework in the Workflow will not start before the currently running 
    Firework in the Workflow completes.

    Parameters:
        -args (objs):       BgwFirework objects
        -deps_dict {dict}:  Specifies the dependency of the BgwInput objects given.
                            If no dependency is given, Fireworks are assumed to be 
                            sequentially dependent.

        -name (str):        Name to be to the Workflow
        -db_upload (bool)   Choose whether or not to upload Successfully completed
                            runs directly to db using bgw_db.yaml defaults file 
                            from HOME Directory

    Example:
        BgwWorkflow(FW1, FW2, FW3, FW4, deps_dict={FW1: [FW2, FW3], 
                        FW2: [FW4], FW3: [FW4]}, name='Example WF')

        This will create a Workflow containing the 4 given BgwFirework objects
        with a Workflow name of 'Example WF' and the given dependencies.
        Dependency Dictionary Explanation:
            In the above example, FW2 and FW3 will not start before FW1 is complete.
            Likewise, FW4 depends on the completion of FW2 and FW3 before starting.
    '''

    def __init__(self, *args, **kwargs):
        '''
        :param args:        (VaspFirework objects) objects to create Workflow from.  
                            No limiton the amount of VaspInputInterface objects to 
                            be given.  Entered as just comma separated objects 
                            passed to class.
        :param deps_dict:   (dict) specifies the dependency of the BgwInput objects 
                            given.  If no dependency is given, Firworks are assumed 
                            to be sequentially dependent.
        :param name         (str) Name given to Workflow
        :param db_upload:   (bool) Set to True for automatic insertion to MongoDB
                            database.  Will read 'bgw_db.yaml' file for DB information.
                            Default: False
        '''
        self.fws = []
        self.name = kwargs.get('name', 'Sequential WF')
        self.deps_dict = kwargs.get('deps_dict', {}) 
        self.db_upload = kwargs.get('db_upload', False)
        self.dependency = {}
        if self.deps_dict:
            for i in self.deps_dict.keys():
                fw_deps = []
                for j in self.deps_dict[i]:
                    ( fw.deps.append(j) if isinstance(j, Firework) \
                            else fw_deps.append(j.Firework) )
                if isinstance(i, Firework): 
                    self.dependency[i] = fw_deps
                else: 
                    self.dependency[i.Firework] = fw_deps 
        self.deps = True if self.dependency else False
        for id, fw_task in enumerate(args):
            self.fws.append(fw_task) if isinstance(fw_task, 
                    Firework) else self.fws.append(fw_task.Firework)
            if not self.deps and id != 0:
                    self.dependency[self.fws[id-1]]=( [fw_task] if isinstance(fw_task,
                                        Firework) else [fw_task.Firework] )
        if self.db_upload:
            db_fw = Firework(BgwDB(config_file='bgw_db.yaml', insert_to_db=True), name="BGW DB Task")
            self.fws.append(db_fw)
            self.dependency[self.fws[id]] = [db_fw]
        self.wf = Workflow(self.fws, self.dependency, name=self.name)

        # Try to establish connection with Launchpad
        try:
            self.LaunchPad=LaunchPad.from_file(os.path.join(
                os.environ["HOME"],".fireworks", "my_launchpad.yaml"))
        except:
            self.LaunchPad = None

    def add_fw(self, fw_task, deps=None):
        self.fws.append(fw_task) if isinstance(fw_task,
                Firework) else self.fws.append(fw_task.Firework)
        if deps:
            for i in deps.keys():
                fw_deps = [j if isinstance(j, Firework) else j.Firework for j in deps[i]]
                
                if not isinstance(i, Firework):
                    i = i.Firework
   
                if i in self.dependency.keys():
                    self.dependency[i].extend(fw_deps)
                else:
                    self.dependency[i]=fw_deps
    
        else:
            id = len(self.fws) - 2
            self.dependency[self.fws[id]]=[fw_task.Firework]
        self.wf=Workflow(self.fws, self.dependency, name=self.name)


    def to_file(self, filename):
        self.wf.to_file(filename)


    def preserve_worker(self):
        self.add_spec('_preserve_fworker', True)


    def add_spec(self, key, val):
        for i in self.fws:
            i.spec[key] = val

    def add_handler(self, handler, **kwargs):
        '''
        Member function to add handler and handler options to
        all Fireworks in the BgwWorkflow.BgwFirework.tasks list.
        '''
        for i in self.fws:
            for j,k in enumerate(i.tasks):
                if instance(k, BgwCustodianTask):
                    cust_handlers = k.get('handlers', []) 
                    cust_params = k.get('handler_params', {})
                    if handler not in cust_handlers:
                        k['handlers'].append(handler)
                        k['handler_params'].update(kwargs)
            i.spec['_tasks'] = [t.to_dict() for t in
                    i.tasks] #re-initialize FW.spec


    def add_wf_to_launchpad(self):
        if self.LaunchPad:
            self.LaunchPad.add_wf(self.wf)
        else:
            print("No connection to LaunchPad. \n"
                    "Use 'to_file(<filename>)' to write a yaml file\n"
                    "to manually add Workflow to LaunchPad later.\n")

