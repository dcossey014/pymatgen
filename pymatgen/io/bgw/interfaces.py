import os, sys
import glob
import shutil
from math import ceil

import subprocess
from monty.serialization import loadfn
from fireworks import Firework, FireTaskBase, FWAction, \
                        explicit_serialize, FileTransferTask, Workflow, LaunchPad
from custodian import Custodian
from pymatgen import Structure
from pymatgen.io.bgw.kgrid import QeMeanFieldGrids, generate_kpath
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

@explicit_serialize
class QeMeanFieldTask(FireTaskBase):
    #TODO: Create to_file()
    '''
    Custodian Task for running Mean Field Calculations with 
    Quantum Espresso.
    ''' 
    required_params = ['structure', 'pseudo_dir', 'kpoints_coarse',
                    'kpoints_fine', 'qshift', 'fftw_grid', 
                    'mpi_cmd', 'pw_cmd', 'pw2bgw_cmd', 'mf_tasks']

    optional_params = ['control', 'electrons', 'system', 'cell', 'ions',
                    'reduce_structure', 'bgw_rev_off', 'log_cart_kpts', 
                    'pw2bgw_input', 'bandstructure_kpoint_path',
                    'num_kpoints_bandstructure', 'config_file']

    def __init__(self, structure, pseudo_dir=None, kpoints_coarse=None,
                kpoints_fine=None, qshift=None, fftw_grid=None, mpi_cmd='', 
                pw_cmd='', pw2bgw_cmd='', mf_tasks=[], control={}, 
                electrons={}, system={}, cell={}, ions={}, 
                reduce_structure=False, bgw_rev_off=False, 
                log_cart_kpts=False, pw2bgw_input={}, 
                bandstructure_kpoint_path=False, num_kpoints_bandstructure=500,
                config_file=None):

        # Mandatory Parameters
        self.structure = structure
        self.pseudo_dir = pseudo_dir
        self.kpoints_co = kpoints_coarse
        self.kpoints_fi = kpoints_fine
        self.qshift = qshift
        self.fftw_grid = fftw_grid
        self.mpi_cmd = mpi_cmd
        self.pw_cmd = pw_cmd
        self.pw2bgw_cmd = pw2bgw_cmd
        self.mf_tasks = mf_tasks

        print("in QEMFT: structure: {}\npseudo: {}\n, kpt_co: {}\nkpt_fi: {}".format(
                self.structure, self.pseudo_dir, self.kpoints_co, self.kpoints_fi))
        # Specific Dictionaries for each calculation
        self.control = control or {}
        self.electrons = electrons or {}
        self.system = system or {}
        self.cell = cell or {}
        self.ions = ions or {}
        self.pw2bgw_input = pw2bgw_input or {}

        # Other Parameters
        self.reduce_structure = reduce_structure or False
        self.bgw_rev_off = bgw_rev_off or False
        self.log_cart_kpts = log_cart_kpts or False
        self.bandstructure_kpoint_path = bandstructure_kpoint_path or  False
        self.num_kpoints_bandstructure = num_kpoints_bandstructure or 500
        self.config_file = config_file if config_file else None

        # Get Primitive Sructure if reducing structure
        if self.reduce_structure:
            self.structure = self.structure.get_primitive_structure()

        # Get Kgrids for QEMF
        if self.bandstructure_kpoint_path:
            self.mf_tasks = ['scf', 'wfn_fi']
            kpath = Generate_Kpath(self.structure, self.num_kpoints_bandstructure)
            kpoints_fi = kpath.create_path()
        
        print("in QEMFT: structure: {}".format(self.structure))
        qemf_kgrids = QeMeanFieldGrids(self.structure, kpoints_coarse=self.kpoints_co,
                            kpoints_fine=self.kpoints_fi, qshift=self.qshift,
                            fftw_grid=self.fftw_grid, bgw_rev_off=self.bgw_rev_off, 
                            log_cart_kpts=self.log_cart_kpts)

        self.kgrids = qemf_kgrids.generate_kgrids()


        # Use QeMFInputs from Espresso Inputs to build inputs.
        self.inputs = QeMFInputs(self.structure, pseudo_dir = self.pseudo_dir,
                        kpoints_coarse=self.kgrids, kpoints_fine=self.kgrids,
                        config_file=self.config_file)

        for i in self.mf_tasks:
            for j in ['control', 'system', 'cell', 'ions', 'electrons']:
                d = self.get(j, {}).get(i, {})
                if d:
                    task = self.inputs.get(i, {}).get(j, {})
                    task.update(d)
                    
    # Run each Espresso Calculation
    def run_task(self, fw_spec):
        self.inputs.to_file()
        
        # Method for running each Espresso Task
        def run_qe_task(self, task):
            os.chdir(task)
            if qe_task == 'scf' :
                self.run_pw(self.pw, qe_task)
            else:
                self.run_pw(self.pw, qe_task, pw2bgwx=self.pw2bgw)
            self.qe_dirs[qe_task] = os.getcwd()
            os.chdir("../")

        for task in self.mf_tasks:
            run_qe_task(task)
        os.chdir("../")
        return FWAction(stored_data=self.output, mod_spec=[
                        {'_set': {'PREV_DIRS': self.prev_dirs}}])

    '''
    def run_task(self, fw_spec):

        # set up directory structure in ./ESPRESSO and create inputs and links
        # and set up some class variables
        self.write_inputs()

        def run_qemf_task(qe_task):
            os.chdir(qe_task)
            if qe_task == 'scf' :
                self.run_pw(self.pw, qe_task)
            else:
                self.run_pw(self.pw, qe_task, pw2bgwx=self.pw2bgw)
            self.qe_dirs[qe_task] = os.getcwd()
            os.chdir("../")

        os.chdir('./ESPRESSO')

        # run all Quantum Espresso mean field tasks
        for qe_task in  self.qe_tasks:
            run_qemf_task(qe_task)
        
        os.chdir("../")
        return FWAction(stored_data=self.output, mod_spec=[{'_set': {'PREV_DIRS': self.prev_dirs}}])


    def write_inputs(self):
        self.structure = self.get('structure')
        self.kpoints = self.get('kpoints')
        self.pseudo_dir = self.get('pseudo_dir')
        self.mpi_cmd = self.get('mpi_cmd', None).split()
        kgrid_offset_type = self.get('krgrid_offset_type', 'Monkhorst-Pack')
        qshift = self.get('qshift', [0, 0, 0.001])
        fftw_grid = self.get('fftw_grid', [0,0,0])
        bgw_rev_off = self.get('bgw_rev_off', False)
        log_cart_kpts = self.get('log_cart_kpts', False)
        self.alternate_kpoints = self.get('bandstructure_kpoint_path', False)
        self.num_kpoints = self.get('num_kpoints', 500)
        self.pw = self.get('pw_cmd', 'pw.x').split()
        self.pw2bgw = self.get('pw2bgw_cmd', 'pw2bgw.x').split()
        self.prev_dirs = {'ESPRESSO': {}}
        self.qe_dirs = self.prev_dirs['ESPRESSO']
        self.output = {}

        #Set PW input Sections
        qe_control = self.get('qe_control', {})
        qe_system = self.get('qe_system', {})
        qe_electrons = self.get('qe_electrons', {})
        qe_pw2bgw = self.get('qe_pw2bgw', {})

        mf_tasks = self.get('mf_tasks', {})
        self.qe_tasks=mf_tasks['qe_tasks']

        # Set Pseudopotential dictionary and gather PP files
        def either(c):
            return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
        self.pseudo = {}
        for i in self.structure.symbol_set:
            pattern="{}/{}_*.UPF".format(self.pseudo_dir,i)
            new_pattern=''.join(either(char) for char in pattern)
            pps=glob.glob(new_pattern)
            self.pseudo[i.encode('ascii', 'ignore')] = pps[-1].split('/')[-1]
            

        # Set up basic kgrid object from which various grids are derived for
        # each mean field task
        qe_kgrids = QeMeanFieldGrids(self.structure, kpoints=self.kpoints,
                            offset_type=kgrid_offset_type, qshift=qshift,
                            fftw_grid=fftw_grid, bgw_rev_off=bgw_rev_off,
                            log_cart_kpts=log_cart_kpts)

        if not os.path.exists('./ESPRESSO'):
            os.mkdir('./ESPRESSO')
        os.chdir('./ESPRESSO')

        # function for setting up an mean field task
        #gk: adjustments should be avoided at runtime. They should be done
        #gk: when the fw spec is created
        def make_qemf_input(qe_task):
            self.dir_setup(qe_task)
            if self.alternate_kpoints and qe_task != 'scf':
                gen_kpoints = generate_kpath(self.structure, 
                                self.num_kpoints)
                qe_task_grid = ['{:5d}\n'.format(len(gen_kpoints))]
                qe_task_grid.extend(
                        [' {:15.10f} {:15.10f} {:15.10f}   1.0\n'.format(
                            k[0], k[1], k[2]) for k in gen_kpoints] )
                qe_task_grid[-1] = qe_task_grid[-1].rstrip()
            else:
                qe_task_grid = qe_kgrids.generate_kgrid(qe_task)
        
            qe_task_control = qe_control.get(qe_task)

            qe_task_system  = qe_system.get(qe_task)
            qe_task_electrons = qe_electrons.get(qe_task)

            # Get number of Bands needed for a good Calculation
            #gk: this logic is not correct. Some cases use number of occupied 
            #gk:  bands: wnfq and wfnq_fi 
            bands = 0
    	    composition = self.structure.composition.as_dict()
            for element in self.pseudo.keys():
                with open(os.path.join("{}/{}".format(
                    self.pseudo_dir,self.pseudo[element])) ) as fin:
                    for line in fin:
                        if "z_valence" in line:
                            l = line.strip().split()
                            z_val = float(l[-1][:-1])
                            bands += (z_val / 2) * composition[element]
                            break

            bands = int(ceil(bands + 4)) if bands * 0.2 < 4 else int(ceil(bands * 1.2))

            # Check if nbnd is set and adjust if necessary
            if 'nbnd' in qe_task_system.keys():
                if qe_task_system['nbnd'] < bands:
                    print("Setting number of bands to: {}\n"
                            "Number of bands was set less "
                            "than the recommended number of bands".format(bands))
                    qe_task_system['nbnd'] = bands
            else:
                qe_task_system['nbnd'] = bands

            #Write Input file for specified QE PWSCF Mean field task
            qe_task_pw = PWInput(self.structure, control=qe_task_control,
                         pseudo=self.pseudo, system=qe_task_system,
                         electrons=qe_task_electrons, kpoints_grid=qe_task_grid, 
                         kpoints_mode='crystal')
            qe_task_pw.write_file('in')

            #Write input file for pw2bgw to convert to wave function to BGW format
            if qe_task != 'scf' and not self.alternate_kpoints: 
                print 'qe_task = ', qe_task
                qe_task_pw2bgw = qe_pw2bgw.get(qe_task)
                #print 'qe_task_pw2bgw = ', qe_task_pw2bgw
                #print 'qe_task_grid = ', qe_task_grid
                qe_task_grid_kpoints=qe_kgrids.grids[qe_task].kpoints
                print 'qe_task_grid_kpoints = ', qe_task_grid_kpoints
                qe_task_grid_offset_type=qe_kgrids.grids[qe_task].offset_type
                print 'qe_task_grid_offset_type = ', qe_task_grid_offset_type
                qe_task_grid_qshift=qe_kgrids.grids[qe_task].qshift
                print 'qe_task_grid_qshift = ', qe_task_grid_qshift
                qe_pw2bgw_input = Pw2BgwInput(self.structure, 
                             pw2bgw_input=qe_task_pw2bgw, 
                             kpoints=qe_task_grid_kpoints,
                             k_offset=qe_task_grid_offset_type,
                             qshift=qe_task_grid_qshift
                             )
                qe_pw2bgw_input.write_input('pp_in')

            os.chdir("../")

        # Do the work setting up the input now for all tasks
        for qe_task in  self.qe_tasks:
            make_qemf_input(qe_task)

        os.chdir("../")


    def dir_setup(self, dir):
        save_dir = '{}.save'.format(self.structure.formula.replace(' ', ''))
        if not os.path.exists(dir):
            os.mkdir(dir)

        os.chdir(dir)
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        if not 'scf' in dir:
            if not os.path.lexists('{}/data-file.xml'.format(save_dir)):
                os.symlink('../../scf/{}/data-file.xml'.format(save_dir),
                        '{}/data-file.xml'.format(save_dir))
            if not os.path.lexists('{}/charge-density.dat'.format(save_dir)):
                os.symlink('../../scf/{}/charge-density.dat'.format(save_dir),
                        '{}/charge-density.dat'.format(save_dir))
    '''

    def run_pw(self, pwx, dir_name, pw2bgwx=None):
        mpi_pw = list(self.mpi_cmd)
        mpi_pw.extend(pwx)
        if pw2bgwx and not self.alternate_kpoints:
            mpi_pw2bgw = list(self.mpi_cmd)
            mpi_pw2bgw.extend(pw2bgwx)
            job = PWJob(mpi_pw, pw2bgw_cmd=mpi_pw2bgw)
        else:
            job = PWJob(mpi_pw)
        handlers = []
        handlers.append(load_class("pymatgen.io.bgw.handlers", 'QuantumEspressoErrorHandler')())
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

