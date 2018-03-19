import os, sys, errno
import glob
import shutil
from math import ceil
from copy import deepcopy as dc

import subprocess
from monty.serialization import loadfn
from fireworks import Firework, FireTaskBase, FWAction, \
                        explicit_serialize, FileTransferTask, Workflow, LaunchPad
from custodian import Custodian
from pymatgen import Structure
from pymatgen.io.espresso.inputs import QeMFInput, QeMFPw2BgwInputs
from pymatgen.io.bgw.kgrid import QeMeanFieldGrids, generate_kpath
from pymatgen.io.bgw.inputs import BgwInput, BgwInputTask
from pymatgen.io.bgw.custodian_jobs import BGWJob, BgwDB, BgwCustodianTask
from pymatgen.io.espresso.custodian_jobs import PWJob

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

        print "gk: in BgwFirework() __init__"
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
        runtype = bgw_task.params['in_file'].split('/')[-1].split('.')[0]
        self.in_file = "OUT.{}".format(runtype[:3])
        self.bgw_cmd = bgw_cmd
        self.ppx = ppx 
        if mpi_cmd:
            self.mpi_cmd = mpi_cmd
        else:
            print("Error: No mpi_cmd given. "
                    "Please input an mpi_cmd for the job.  \n"
                    "Usage: bgw_fw = BgwFirework(<bgw_task>, mpi_cmd=<mpi_cmd>, "
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
                fout=self.in_file, mpi_cmd=self.mpi_cmd, 
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

@explicit_serialize
class QeMeanFieldTask(FireTaskBase):
    #TODO: Create to_file()
    #TODO: Create as_dict() and from_dict() functions
    '''
    Custodian Task for running Mean Field Calculations with 
    Quantum Espresso and postprocessed for use with BerkeleyGW.
    ''' 
    required_params = ['structure', 'kpoints_coarse', 
                    'kpoints_fine', 'qshift', 'fftw_grid', 
                    'mpi_cmd', 'pw_cmd', 'pw2bgw_cmd', 'mf_tasks']

    optional_params = ['control', 'electrons', 'system', 'cell', 'ions',
                    'reduce_structure', 'bgw_rev_off', 'log_cart_kpts', 
                    'pw2bgw_input', 'bandstructure_kpoint_path',
                    'num_kpoints_bandstructure', 'cmplx_real', 
                    'config_file','pseudo_dir']
    
    def __init__(self, *args, **kwargs):
        print "gk: in __init__ QeMeanFieldTask in interfaces in bgw dir"
        super(QeMeanFieldTask, self).__init__(*args, **kwargs)
        print "gk: leaving __init__ QeMeanFieldTask in interfaces in bgw dir"

    def build_inputs(self, dry_run=False):
        # Mandatory Parameters
        self.structure = self.get('structure')
        self.pseudo_dir = self.setdefault('pseudo_dir',{})
        self.kpoints_co = self.get('kpoints_coarse')
        self.kpoints_fi = self.get('kpoints_fine')
        self.qshift = self.get('qshift')
        self.fftw_grid = self.get('fftw_grid')
        self.mpi_cmd = self.get('mpi_cmd').split()
        self.pw_cmd = self.get('pw_cmd').split()
        self.pw2bgw_cmd = self.get('pw2bgw_cmd').split()
        self.mf_tasks = self.get('mf_tasks')

        # Other Parameters
        self.reduce_structure = self.get('reduce_structure', False)
        self.bgw_rev_off = self.get('bgw_rev_off', False)
        self.log_cart_kpts = self.get('log_cart_kpts', False)
        self.bandstructure_kpoint_path = self.get('bandstructure_kpoint_path',  False)
        self.num_kpoints_bandstructure = self.get('num_kpoints_bandstructure', 500)
        self.cmplx_real = self.get('cmplx_real', 'real')
        self.config_file = self.get('config_file', None)

        # Import Configuration Dictionary if given a configuration File.
        if self.config_file and os.path.isabs(self.config_file):
            self.__dict__['config_dict'] = loadfn(self.config_file)
        elif self.config_file:
            fpath = os.path.join(os.environ['HOME'], self.config_file)
            self.__dict__['config_dict'] = loadfn(fpath)
        elif os.path.exists(os.path.join(os.environ['HOME'],
                    'bgw_interface_defaults.yaml')):
            fpath = os.path.join(os.environ['HOME'], 'bgw_interface_defaults.yaml')
            self.__dict__['config_dict'] = loadfn(fpath)
        else:
            self.__dict__['config_dict'] = {}


        # Specific Dictionaries for each calculation
        prefix = self.get('control', {}).get('scf', {}).get(
                            'prefix', '')
        prefix = ( self.structure.composition.reduced_formula 
                        if not prefix else prefix )
        espresso_dict = self.config_dict.get('espresso', {} )
        self.pw2bgw_input_dict = self.config_dict.get('pw2bgw', {} )
        wfn_file = "wfn.{}".format(self.cmplx_real.lower())
        rhog_file = "rho.{}".format(self.cmplx_real.lower())
        cmplx_real_num = 1 if 'real' in self.cmplx_real.lower() else 2

        (self.__dict__['task_control'], 
            self.__dict__['task_system'],
            self.__dict__['task_electrons'],
            self.__dict__['task_ions'],
            self.__dict__['task_cell'] ) = {}, {}, {}, {}, {}

        if self.bandstructure_kpoint_path:
            self.mf_tasks.append('wfn_bs')

        for i in self.mf_tasks:
            # Take in Defaults from User YAML

            self.task_control[i] = dc(espresso_dict.get('control', {}) )
            self.task_control[i].update({'verbosity': 'high', 'prefix': prefix})

            self.task_system[i] = dc(espresso_dict.get('system', {}) )
            self.task_electrons[i] = dc(espresso_dict.get('electrons', {}) )
            self.task_ions[i] = dc(espresso_dict.get('ions', {}) )
            self.task_cell[i] = dc(espresso_dict.get('cell', {}) )
            if 'scf' not in i:
                pw2bgw_task_dict = self.pw2bgw_input_dict.get(i, {})
                pw2bgw_task_dict.update({'prefix': prefix, 'wfng_file': wfn_file,
                                            'real_or_complex': cmplx_real_num} 
                                        )

            # Update with User Inputed Values
            self.task_control[i].update(self.get('control', {}).get(i, {}) )
            self.task_system[i].update(self.get('system', {}).get(i, {}) )
            self.task_electrons[i].update(self.get('electrons', {}).get(i, {}) )
            self.task_ions[i].update(self.get('ions', {}).get(i, {}) )
            self.task_cell[i].update(self.get('cell', {}).get(i, {}) )

        print("System: {}".format(self.task_system))
        print("Control: {}".format(self.task_control))
        print("electrons: {}".format(self.task_electrons))

        # Set VXC_Diag_Nmax and Update Pw2Bgw_Input dictionary with User values
        print("pw2bgw_dict: {}".format(self.pw2bgw_input_dict))
        wfn_co_pw2bgw = self.get('pw2bgw_input_dict', {}).get('wfn_co', {})
        print("wfn_co_pw2bgw: {}".format(wfn_co_pw2bgw))
        wfn_co_pw2bgw['vxc_diag_nmax'] = self.task_system['wfn_co']['nbnd']
        wfn_co_pw2bgw['rhog_file'] = rhog_file
        self.pw2bgw_input_dict.update(self.get('pw2bgw_input', {}) )

        # Get Primitive Sructure if reducing structure
        if self.reduce_structure:
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            finder = SpacegroupAnalyzer(self.structure)
            self.__dict__['structure'] = finder.get_primitive_standard_structure()

        qemf_kgrids = QeMeanFieldGrids(self.structure, kpoints_coarse=self.kpoints_co,
                            kpoints_fine=self.kpoints_fi, qshift=self.qshift,
                            fftw_grid=self.fftw_grid, bgw_rev_off=self.bgw_rev_off, 
                            log_cart_kpts=self.log_cart_kpts)

        self.__dict__['kgrids'] = qemf_kgrids.generate_kgrids()

        # Get Kgrids for QEMF Bandstructure Plots if given bandstructure_kpoint_path
        if self.bandstructure_kpoint_path:
            kpath = Generate_Kpath(self.structure, self.num_kpoints_bandstructure)
            self.__dict__['kgrids']['wfn_bs'] = kpath.create_path()
        

        # Use QeMFInputs from Espresso Inputs to build inputs.
        self.__dict__['inputs'] = QeMFInput(self.structure, pseudo_dir = self.pseudo_dir,
                        kpoints_coarse=self.kgrids, kpoints_fine=self.kgrids,
                        control=self.task_control, system=self.task_system, 
                        ions=self.task_ions, electrons=self.task_electrons, 
                        cell=self.task_cell, tasks=self.mf_tasks,
                        reduce_structure=self.reduce_structure) 

        # Use QeMFPw2BgwInputs to build PostProcessing Inputs.
        self.__dict__['pw2bgw_inputs'] = QeMFPw2BgwInputs(self.structure,
                        pw2bgw_input=self.pw2bgw_input_dict, mf_tasks=self.mf_tasks)

        # Remove ESPRESSO directory from __init__.qemf_kgrids if dry_run=True for debugging
        if dry_run:
            shutil.rmtree('./ESPRESSO')

    def write_inputs(self):
        # Write Input Files
        self.inputs.to_file()
        self.pw2bgw_inputs.to_file()
        
        # Symbolically link data file and charge density file from SCF
        #   to all other tasks.
        def mkdir(dir):
            if not os.path.exists(dir):
                os.mkdir(dir)

        def force_link(file1, file2):
            try:
                os.symlink(file1,file2)
            except OSError, e:
                if e.errno == errno.EEXIST:
                    os.remove(file2)
                    os.symlink(file1,file2)

        prefix  = self.inputs.scf.control['prefix']
        save_dir = "{}.save".format(prefix)
        scf_dir = "ESPRESSO/scf/{}".format(save_dir)
        scf_chrg = os.path.join('../../../', scf_dir, 'charge-density.dat')
        scf_dat = os.path.join('../../../', scf_dir, 'data-file.xml')

        mkdir(scf_dir)
        for i in self.mf_tasks:
            if 'scf' not in i.lower():
                link_dir = os.path.join('ESPRESSO', i, save_dir)
                mkdir(link_dir)
                link_chrg = os.path.join(link_dir, 'charge-density.dat')
                link_dat = os.path.join(link_dir, 'data-file.xml')

                force_link(scf_dat, link_dat)
                force_link(scf_chrg, link_chrg)
                    
    # Run each Espresso Calculation
    def run_task(self, fw_spec):
        # Write Input files and make Symbolic links
        self.build_inputs()
        self.write_inputs()
        print("mpi_cmd: {}\npw_cmd: {}\npw2bgw_cmd: {}".format(
                self.mpi_cmd, self.pw_cmd, self.pw2bgw_cmd) )
        #sys.exit(0)
        
        # Method for running each Espresso Task
        def run_qe_task(task):
            os.chdir(task)
            if task == 'scf' :
                self.run_pw(self.pw_cmd, task)
            else:
                self.run_pw(self.pw_cmd, task, pw2bgwx=self.pw2bgw_cmd)
            self.prev_dirs['ESPRESSO'][task] = os.getcwd()
            os.chdir("../")

        # Change to Top Level Directory and begin running pw.x
        ( self.__dict__['output'],
            self.__dict__['prev_dirs'] ) = {}, {'ESPRESSO': {} }
        os.chdir("./ESPRESSO")
        for task in self.mf_tasks:
            run_qe_task(task)
        os.chdir("../")
        return FWAction(stored_data=self.output, mod_spec=[
                        {'_set': {'PREV_DIRS': self.prev_dirs}}])

    def run_pw(self, pwx, dir_name, pw2bgwx=None):
        mpi_pw = list(self.mpi_cmd)
        mpi_pw.extend(pwx)
        if pw2bgwx and not self.bandstructure_kpoint_path:
            mpi_pw2bgw = list(self.mpi_cmd)
            mpi_pw2bgw.extend(pw2bgwx)
            job = PWJob(mpi_pw, pw2bgw_cmd=mpi_pw2bgw)
        else:
            job = PWJob(mpi_pw)
        handlers = []
        handlers.append(load_class("pymatgen.io.espresso.handlers", 'QuantumEspressoErrorHandler')())
        c = Custodian(handlers=handlers, validators=[], jobs=[job])
        self.output[dir_name] = c.run()

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

