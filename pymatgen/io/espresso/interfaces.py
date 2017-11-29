import os, sys, glob, shutil, subprocess
from math import ceil

from monty.serialization import loadfn
from fireworks import Firework, FireTaskBase, FWAction, \
                        explicit_serialize, FileTransferTask, Workflow, LaunchPad
from custodian import Custodian
from pymatgen import Structure
from pymatgen.io.bgw.kgrid import QeMeanFieldGrids, generate_kpath
from pymatgen.io.espresso.inputs import PWInput, Pw2BgwInput
from pymatgen.io.espresso.custodian_jobs import PWJob

def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)


@explicit_serialize
class QeMeanFieldTask(FireTaskBase):
    #TODO: Create to_file()
    '''
    Custodian Task for running Mean Field Calculations with 
    Quantum Espresso.
    ''' 
    required_params = ['structure', 'pseudo_dir', 'kpoints_coarse',
                    'kpoints_fine', 'mpi_cmd', 'pw_cmd', 'pw2bgw_cmd',
                    'mf_tasks']
    optional_params = ['reduce_structure', 'bgw_rev_off', 'log_cart_kpts', 
                    'pw_cmd', 'pw2bgw_cmd', 'bandstructure_kpoint_path', 
                    'num_kpoints']

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


