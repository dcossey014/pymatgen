__author__ = "David Cossey"
__copyright__ = "n/a"
__version__ = "0.1alpha"
__maintainer__ = "David Cossey"
__email__ = "dcossey014@gmail.com"
__date__ = "2016-04-28"

import six, glob, sys, errno
import os, abc
import mmap, fnmatch, re, subprocess
import copy as cp
from bisect import bisect_left

from pymatgen.io.bgw.kgrid import Kgrid, generate_kpath

from monty.io import zopen
from monty.dev import deprecated
from monty.json import MSONable
from monty.serialization import loadfn
from collections import defaultdict, OrderedDict
from fireworks import FireTaskBase, explicit_serialize, FWAction


@explicit_serialize
class BgwInputTask(FireTaskBase):
    '''
    Fireworks Task for writing BerkeleyGW Input Files.

    Required Parameters:
        input_set_params (dict):Dictionary of input set parameters for the BerkeleyGW
                                job to be run
        in_file    (str)   :   Filename to write Input file to

    Optional Parameters:
        kpoints (array)     :   List of Kpoints used in Kgrid.x to create WFN.
                                This is required for Epsilon and Sigma calculations.
        qshift  (array)     :   Qshift used in Kgrid.x for WFN Qshifted grids.
                                This is required for Epsilon calculations.
        mat_type (str)      :   Options: 'metal' or 'semiconductor'.  This is required
                                for Epsilon calculations.

    '''
    
    required_params = ['input_set_params', 'in_file']
    optional_params = ['kpoints', 'qshift', 'mat_type', 'cmplx_real', 'qemf_dir', 'config_file', 
        'reduce_structure', 'input_set_params', 'in_file']

    def __init__(self, params): 
        print "\n\nin BgwInputTask __init__"

        self.kpoints = params.get('kpoints', None)
        self.qshift = params.get('qshift', None)
        self.isp = params.get('input_set_params')
        self.run_type = params.get('run_type', None)
        self.mat_type = params.get('mat_type', 'semiconductor')
        self.filename = params.get('in_file')
        self.cmplx_real = params.get('cmplx_real')
        self.qemf_dir = params.get('qemf_dir')
        self.kps = params.get('kps', None)
        #gk: why is occupied_bands set here?
        #self.occupied_bands = params.get('occupied_bands', 0)
        self.config_file = params.get('config_file')
        self.reduce_structure = params.get('reduce_structure', False)

        params = {
                'qshift': self.qshift, 
                'input_set_params': self.isp, 'run_type': self.run_type,
                'mat_type':self.mat_type, 'in_file': self.filename,
                'kps': self.kps, 'cmplx_real': self.cmplx_real,'qemf_dir': self.qemf_dir}

        self.update(params)

    def convert2primitive(self, s):
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        finder = SpacegroupAnalyzer(s)
        return finder.get_primitive_standard_structure()

    def __str__(self):
        out = []
        def write_kpoints():
            if not self.kps:
                kp_dir = self.prev_dirs['ESPRESSO']['wfn_co']
                q_shift_dir = self.prev_dirs['ESPRESSO']['wfnq']
                wfn_out = os.path.join(kp_dir, 'wfn_co.out')

                with open(os.path.join(q_shift_dir, 'wfnq.in'), 'r') as fin:
                    lines = fin.readlines()
                
                qshift = [float(i) for i in lines[2].strip().split()]
                self.kps = Kgrid.from_file(wfn_out)

            qtype = 2 if self.mat_type == 'metal' else 1

            if 'epsilon' in self.run_type.lower():
                self.kps = [[j[0], j[1], j[2], 1.0, 0] for j in self.kps]
                self.kps[0] = [qshift[0], qshift[1], qshift[2], 1.0, qtype]

                out.extend(["  {:.6f}  {:.6f}  {:.6f} {:>4.1f}  {:2d}".format(
                                j[0], j[1], j[2], j[3], j[4]
                                ) for j in self.kps ]
                            )
            else:
                out.extend(["  {:.6f}  {:.6f}  {:.6f} {:>3.1f}".format(
                                j[0], j[1], j[2], 1.0
                                ) for j in self.kps ]
                            )

            out.append('end')

        #def write_band_occ():
        #    print "gk in BgwInput.__str__.write_band_occ(), self.occupied_bands=",self.occupied_bands
        #    out.append("band_occupation {}*1 {}*0".format(
        #            int(self.occupied_bands), 
        #            int(self.isp['number_bands']) - int(self.occupied_bands)))

        for k1 in self.isp.keys():
            out.append("{}  {}".format(k1, self.isp[k1]))
        print "gk: first out = ", out

        if 'epsilon' in self.run_type:
            #if self.occupied_bands != None:
            #    write_band_occ()
            out.append('begin qpoints')
            write_kpoints()
        elif 'sigma' in self.run_type:
            #if self.occupied_bands != None:
            #    write_band_occ()
            out.append('begin kpoints')
            write_kpoints()

        return('\n'.join(out))
    
    def write_file(self, filename):
        '''
        Write input file for BGW.
        
        Args:
            filename(str):  The string filename to output to.
        '''
        if not self.run_type:
            self.run_type = filename
        #gk: probably remove below, because now we are checking at FW creation
        print "gk: in BgwInputTask.write_file, filename=",filename
        print "gk: in BgwInputTask.write_file before check_params"
        self.check_params()
        print "gk: in BgwInputTask.write_file after check_params"
        with open(filename, 'w') as f:
            f.write(self.__str__()+'\n')

    def check_params(self):
        print "gk: in BgwInputTask.check_params"
        #gk: probably none of these are required in BGW v 1.2
        eps_params = ['epsilon_cutoff', 'number_bands']
        #gk: probably none of these are required
        #sig_params = ['screened_coulomb_cutoff', 'bare_coulomb_cutoff',
        #            'number_bands', 'band_index_min', 'band_index_max']
        sig_params = ['number_bands', 'band_index_min', 'band_index_max']
        krn_params = ['number_val_bands', 'number_cond_bands', 
                    'screened_coulomb_cutoff']
        abs_params = ['number_val_bands_coarse', 'number_val_bands_fine',
                    'number_cond_bands_coarse', 'number_cond_bands_fine',
                    'energy_resolution']
        sig2wan_params = ['sigma_out', 'wannier_input', 'spin_component',
                    'eqp_num', 'num_bands']

        if 'epsilon' in self.run_type:
            params = list(eps_params)
        elif 'sigma' in self.run_type:
            params = list(sig_params)
            self.isp['band_index_min'] = self.isp.get('band_index_min',1)
            bim=self.isp.get('band_index_max',self.isp.get('number_bands', 0))
            bm=self.isp.get('number_bands')
            bim=bm if bim>bm else bim
            self.isp['band_index_max']=bim
        elif 'sig2wan' in self.run_type:
            params = list(sig2wan_params)
        elif 'abs' in self.run_type:
            params = list(abs_params)
        elif 'kernel' in self.run_type:
            params = list(krn_params)
        else:
            raise NameError('Unrecognized Calculation Type: run_type = %s' 
                        %self.run_type)
            return
        
        for i in params:
            if i not in self.isp.keys():
                raise KeyError('Required Key: "{}"'.format(i) + 
                        ' was not found for given run_type.  ' +
                        'Using run_type: "{}"'.format(self.run_type))
                return

    def dep_setup(self, fout):

        print "gk: in dep_setup"

        def force_link(file1,file2):
            try:
                os.symlink(file1,file2)
            except OSError, e:
                if e.errno == errno.EEXIST:
                    os.remove(file2)
                    os.symlink(file1,file2)

        def extra_links(set):
            if 'eps' in set:
                dir = self.bgw_dirs['epsilon']
                #gk: possibly use logic for non HDF5 case
                #files = ['eps0mat', 'epsmat']
                files = ['eps0mat.h5', 'epsmat.h5']
            elif 'bse' in set:
                dir = self.bgw_dirs['kernel']
                #gk: possibly use logic for non HDF5 case
                #files = ['bsedmat', 'bsexmat']
                files = ['bsemat.h5']
            for f in files:
                force_link(os.path.join(dir, f),
                        os.path.join('.', f))

        print "self.cmplx_real = ", self.cmplx_real
        wfn_file='wfn.'+self.cmplx_real
        wfn_co_file='wfn_co.'+self.cmplx_real

        if 'epsilon' in fout:
            print "Setting up Epsilon dependencies and checking input values."
            wfn = self.qe_dirs['wfn']
            wfnq = self.qe_dirs['wfnq']
            force_link(os.path.join(wfn, wfn_file),
                    './WFN')
            force_link(os.path.join(wfnq, wfn_file),
                    './WFNq')

            self.check_degeneracy(['./WFN'])

        if 'sigma' in fout:
            print "Setting up Sigma dependencies and checking input values."
            wfn_co = self.qe_dirs['wfn_co']
            force_link(os.path.join(wfn_co, 'vxc.dat'),
                    './vxc.dat')
            rho_file='rho.'+self.cmplx_real
            force_link(os.path.join(wfn_co, rho_file),
                    './RHO')
            force_link(os.path.join(wfn_co, wfn_file), 
                    './WFN_inner')
            extra_links('eps')
            self.check_degeneracy(['./WFN_inner'])

        if 'kernel' in fout:
            wfn_co = self.qe_dirs['wfn_co']
            force_link(os.path.join(wfn_co, wfn_file),
                    './WFN_co')
            extra_links('eps')
            self.check_degeneracy(['./WFN_co'])

        if 'absorp' in fout:
            wfn_co = self.qe_dirs['wfn_co']
            wfn_fi = self.qe_dirs['wfn_fi']
            wfnq_fi = self.qe_dirs['wfnq_fi']
            sig_dir = self.bgw_dirs['sigma']
            force_link(os.path.join(wfn_co, wfn_file),
                    './WFN_co')
            force_link(os.path.join(wfn_fi, wfn_file),
                    './WFN_fi')
            force_link(os.path.join(wfnq_fi, wfn_file),
                    './WFNq_fi')
            force_link(os.path.join(sig_dir, 'eqp_co.dat'),
                    './eqp_co.dat')
            extra_links('eps')
            extra_links('bse')
            self.check_degeneracy(['./WFN_fi'])
            self.check_degeneracy(['./WFNq_fi'])
            self.check_degeneracy(['./WFN_co'])

    def check_degeneracy(self, wfn_file, degen_exec='degeneracy_check.x'):
        config_file=self.config_file
        if config_file:
            config_dict = loadfn(config_file)
            de_check_exec = [os.path.join(config_dict['BGW_DIR'], 'degeneracy_check.x')]
        elif os.path.exists(os.path.join(os.environ['HOME'],
                                        'bgw_interface_defaults.yaml')):
            config_dict = loadfn(os.path.join(os.environ['HOME'],
                                        'bgw_interface_defaults.yaml'))
            de_check_exec = [os.path.join(config_dict['BGW_DIR'], 'degeneracy_check.x')]
        else:
            config_dict = {}
            de_check_exec = [degen_exec]

        print "Checking degeneracy using", de_check_exec
        print "Checking wave function file", wfn_file

        de_check_exec.extend(wfn_file)
        p_wfn=subprocess.Popen(de_check_exec,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        p_wfn.wait()

        #gk: check status here and if fails send a message

        allowed_bands, allowed_vbands, allowed_cbands = [], [], []
        lines = p_wfn.stdout.readlines()
        alist = ['buff']
        for line in lines:
            l = line
            print(l)
            if "epsilon" in l:
                alist = allowed_bands
            if "valence" in l:
                alist = allowed_vbands
            if "conduction" in l:
                alist = allowed_cbands

            try:
                if "Note" not in l and l != '\n':
                    alist.append(int(l))
                elif "Note" in l:
                    l = l.split()
                    alist = ['buff']
            except:
                pass
            
        print "    allowed bands=",allowed_bands
        print "    allowed valence bands=",allowed_vbands
        print "    allowed conduction bands=",allowed_cbands

        def match_bands(user_spec,allowed):
            # based on the user specified band number and allowed bands from check_degneracy,
            # match_bands selects the closest degeneracy allowed band number. If there is a 
            # tie, then the larger of the two allowed band numbers will be selected.
            npos = bisect_left(allowed, user_spec)
            if npos == 0:
                return allowed[0]
            if npos == len(allowed):
                return allowed[-1]
            before = allowed[npos -1]
            after = allowed[npos]
            if after - user_spec > user_spec - before:
                return before
            else:
                return after 

        if 'epsilon' in self.filename:
            user_nbands=int(self.isp.get('number_bands', 0))
            print "For Epsilon user specified number_bands =", user_nbands
            if user_nbands < allowed_vbands[-1]:
                print "Number of bands cannot be less than Number of Valence Bands"
                user_nbands = allowed_bands[-1]
                print "Resetting number_bands to, ",user_nbands
            print "Checking for allowed bands closest to number_bands."
            self.allowed_nbands=match_bands(user_nbands,allowed_bands)

            if self.allowed_nbands == user_nbands:
                print "User specified epsilon number_bands is degeneracy allowed"
                self.isp['number_bands'] = user_nbands
            else:
                print "Epsilon number_bands is reset to to degeneracy allowed value,",self.allowed_nbands
                self.isp['number_bands']=self.allowed_nbands

        if 'sigma' in self.filename:
            user_nbands=int(self.isp.get('number_bands', 0))
            print "For Sigma, user specified number_bands =", user_nbands
            if user_nbands < allowed_vbands[-1]:
                print "Number of bands cannot be less than Number of Valence Bands"
                user_nbands = allowed_bands[-1]
                print "Resetting number_bands to, ",user_nbands
            print "Checking for allowed bands closest to number_bands."
            self.allowed_nbands=match_bands(user_nbands,allowed_bands)

            if self.allowed_nbands == user_nbands:
                print "Sigma user specified bands is degeneracy allowed"
                self.isp['number_bands'] = user_nbands
            else:
                print "Sigma number_bands is reset to to degeneracy allowed value,",self.allowed_nbands
                self.isp['number_bands']=self.allowed_nbands

            max_nbands=int(self.isp.get('band_index_max', 0))
            print "For Sigma, user specified band_index_max =", max_nbands
            if max_nbands < allowed_vbands[-1]:
                print "Number of bands cannot be less than Number of Valence Bands"
                exit()
            print "Checking for allowed bands closest to band_index_max."
            self.allowed_nbands=match_bands(max_nbands,allowed_bands)

            if self.allowed_nbands == max_nbands:
                print "Sigma user specified bands is degeneracy allowed"
                self.isp['band_index_max'] = max_nbands
            else:
                print "Sigma band_index_max is reset to to degeneracy allowed value,",self.allowed_nbands
                self.isp['band_index_max']=self.allowed_nbands
                
        if 'kernel' in self.filename:
            user_vband = self.isp.get('number_val_bands', None)
            print("user_vband: {}".format(user_vband))
            self.isp['number_val_bands'] = (allowed_vbands[-1] if not \
                    user_vband else match_bands(user_vband, allowed_vbands) )
            print("self.isp['number_val_bands']: {}".format(self.isp['number_val_bands']))

            user_cband = self.isp.get('number_cond_bands', None)
            self.isp['number_cond_bands'] = (allowed_cbands[-1] if not \
                    user_cband else match_bands(user_cband, allowed_cbands) )

        if 'absorption' in self.filename:
            if 'co' in wfn_file[0]:
                user_vband = self.isp.get('number_val_bands_coarse', None)
                self.isp['number_val_bands_coarse'] = ( allowed_vbands[-1] if not \
                        user_vband else match_bands(user_vband, allowed_vbands))

                user_cband = self.isp.get('number_cond_bands_coarse', None)
                self.isp['number_cond_bands_coarse'] = ( allowed_cbands[-1] if not \
                        user_cband else match_bands(user_cband, allowed_cbands))

            if 'fi' in wfn_file[0] and not 'q' in wfn_file[0]:
                user_cband = self.isp.get('number_cond_bands_fine', None)
                self.isp['number_cond_bands_fine'] = allowed_cbands[-1] if not \
                        user_cband else match_bands(user_cband, allowed_cbands)

            if 'q_fi' in wfn_file[0]:
                user_vband = self.isp.get('number_val_bands_fine', None)
                self.isp['number_val_bands_fine'] = ( allowed_vbands[-1] if not \
                        user_vband else match_bands(user_vband, allowed_vbands))

        print("finished with dep_setup")
        
        
    def dry_run(self, qemf_dir=None):
        self.qemf_dir = qemf_dir if qemf_dir else params.get('qemf_dir')
        print "gk: in BgwInputTask.dry_run(), self.qemf_dir=",self.qemf_dir
        if self.qemf_dir == None:
            print "Can not perform a dry_run without specifying a qemf_dir"
            exit()

        # hard-wire kludge for dry run. user expected to supply this
        self.qe_dirs= {
            "wfnq": self.qemf_dir+"/wfnq",
            "wfn": self.qemf_dir+"/wfn",
            "wfn_fi": self.qemf_dir+"/wfn_fi",
            "scf": self.qemf_dir+"/scf",
            "wfnq_fi": self.qemf_dir+"/wfnq_fi",
            "wfn_co": self.qemf_dir+"/wfn_co"
        }
        self.prev_dirs={ 'ESPRESSO': self.qe_dirs }
        print 'qe_dirs=',self.qe_dirs

        # hard-wire kludge: expecting files to be in current directory
        # hypthetically
        self.bgw_dirs = { "epsilon": "./", "sigma": "./", "kernel": "./" }

        self.dep_setup(self.filename)
        self.write_file(self.filename)
                    


    def run_task(self, fw_spec):
        self.prev_dirs = fw_spec.get('PREV_DIRS', {'ESPRESSO': {}, 'BGW': {}})
        self.qe_dirs = self.prev_dirs.get('ESPRESSO', {})
        if self.qemf_dir:
            self.qe_dirs.update({
                "wfnq": self.qemf_dir+"/wfnq",
                "wfn": self.qemf_dir+"/wfn",
                "wfn_fi": self.qemf_dir+"/wfn_fi",
                "scf": self.qemf_dir+"/scf",
                "wfnq_fi": self.qemf_dir+"/wfnq_fi",
                "wfn_co": self.qemf_dir+"/wfn_co"
                })
            print 'qe_dirs=',self.qe_dirs
            print("prev_dirs: {}".format(self.prev_dirs))

        self.bgw_dirs = self.prev_dirs.get('BGW', {})

        self.dep_setup(self.filename)
        self.write_file(self.filename)
        return FWAction(mod_spec=[{'_set': {
            'PREV_DIRS': self.prev_dirs}}])
                    

class BgwInput(BgwInputTask):
    '''
    Facade Interface for using the BgwInputTask Fireworks Task.

    Parameters:
        filename (str):         Filename to write Input file to.
        kpoints (array):        List of Kpoints used in Kgrid.x to create WFN.
                                This is required for Epsilon and Sigma calculations.
        qshift  (array):        Qshift used in Kgrid.x for WFN Qshifted grids.
                                This is required for Epsilon calculations.
        mat_type (str):         Options: 'metal' or 'semiconductor'.  This is
                                required for Epsilon calculations.
                                Default: semiconductor
    '''

    def __init__(self, isp={}, cmplx_real='cmplx',
                kpoints=None, qshift=None, mat_type='semiconductor',
                kps=None, filename=None,
                qemf_dir=None, reduce_structure=False):

        print "gk: in BgwInput __init__"

        self.__dict__['isp'] = isp if isp else {}
        self.__dict__['run_type'] = filename.split('/')[-1].split('.')[0]
        self.__dict__['params'] = { 
                'kpoints': kpoints, 'qshift': qshift, 'kps': kps,
                'input_set_params': self.isp, 'mat_type': mat_type,
                'in_file': filename, 'run_type': self.run_type,
                'cmplx_real': cmplx_real,
                'qemf_dir': qemf_dir, 'reduce_structure': reduce_structure}

        # this needs to call only with params

        super(BgwInput, self).__init__(self.params)
        print "gk: leaving BgwInput __init__"


    def __setitem__(self, key, val):
        self.proc_key_val(key.strip(), val.strip()) if isinstance(
                val, six.string_types) else self.proc_key_val(key.strip(), val)

    def __setattr__(self, key, val):
        self.proc_key_val(key.strip(), val.strip()) if isinstance(
                val, six.string_types) else self.proc_key_val(key.strip(), val)

    def proc_key_val(self, key, val):
        isp_dict = self.params.get('input_set_params', {})

        valid_params = {
            'epsilon': {
                'int': ['number_bands', 'number_partial_occup', 
                    'number_qpoints', 'frequency_dependence', 
                    'full_chi_cov_log', 'number_valence_pools'],
                'str': ['band_occupation'],
                'bool': ['fermi_level_absolute', 'fermi_level_relative',
                    'eqp_corrections', 'write_vcoul', 'skip_epsilon', 
                    'comm_mpi', 'comm_disk', 'cell_box_truncation',
                    'cell_wire_truncation', 'cell_slab_truncation', 
                    'spherical_truncation', 'gcomm_matrix', 'gcomm_elements',
                    'skip_chi', 'fullbz_replace', 'fullbz_write'
                    'degeneracy_check_override', 'gcomm_matrix', 
                    'gcomm_elements', 'comm_mpi', 'comm_disk'],
                'float': ['epsilon_cutoff', 'fermi_level', 
                    'coulomb_truncation_radius', 'init_frequency', 
                    'delta_frequency', 'delta_frequency_step', 
                    'frequency_low_cutoff', 'frequency_high_cutoff', 
                    'broadening', 'evs', 'ev0', 'evdel', 'ecs', 'ec0', 'ecdel', 
                    'cvfit', 'ecdel_outer', 'cvfit_outer'],
                'exclusive': [['fermi_level_relative','fermi_level_absolute'], 
                    ['ecdel', 'cvfit'], ['ecdel_outer', 'cvfit_outer'],
                    ['comm_mpi', 'comm_disk'], 
                    ['gcomm_matrix', 'gcomm_elements'],
                    ['cell_box_truncation', 'cell_wire_truncation',
                    'cell_slab_truncation', 'spherical_truncation']]
                },
            'sigma': {
                'int': ['number_bands', 'number_partial_occup', 
                    'number_diag', 'begin diag', 'band_index_min', 
                    'band_index_max', 'number_offdiag', 'sigma_matrix', 
                    'spin_index_min', 'spin_index_max', 'frequency_dependence', 
                    'number_frequency_eval', 'interpolation', 
                    'finite_difference_form', 'number_qpoints', 'qgrid',
                    'ggpsum', 'full_ch_conv_log', 'number_kpoints', 
                    'number_sigma_pools', 'exact_static_ch'], 
                'float': ['screened_coulomb_cutoff', 'bare_coulomb_cutoff',
                    'fermi_level', 'coulomb_truncation_radius', 
                    'cell_average_cutoff', 'init_frequency_eval', 
                    'delta_frequency_eval', 'finite_difference_spacing', 
                    'bare_exchange_fraction', 'gpp_broadening', 'gpp_sexcutoff', 
                    'evs', 'ev0', 'evdel', 'ecs', 'ec0', 'ecdel', 'cvfit', 
                    'evs_outer', 'ev0_outer', 'evdel_outer', 'ecs_outer', 
                    'ec0_outer', 'ecdel_outer', 'cvfit_outer', 'avgpot', 
                    'avgpot_outer', 'tol_degeneracy'],
                'str': ['band_occupation'],
                'bool': ['fermi_level_absolute', 'fermi_level_relative', 
                    'screening_metal', 'screening_graphene', 
                    'screening_semiconductor', 'cell_box_truncation', 
                    'cell_wire_truncation', 'cell_slab_truncation', 
                    'spherical_truncation', 'use_epsilon_remainder',
                    'use_xdat', 'dont_use_vxcdat', 'eqp_corrections',
                    'eqp_outer_corrections', 'write_vcoul', 'comm_mpi', 
                    'comm_disk', 'fullbz_replace', 'fullbz_write', 
                    'degeneracy_check_override', 'use_symmetries_q_grid', 
                    'no_symmetries_q_grid', 'die_outside_sphere', 
                    'ignore_outside_sphere', 'no_symmetries_offdiagonals'],
                'exclusive': [['fermi_level_absolute', 'fermi_level_relative'],
                    ['comm_mpi', 'comm_disk'], ['screening_metal', 
                    'screening_graphene', 'screening_semiconductor'], 
                    ['cell_box_truncation', 'cell_wire_truncation',
                        'cell_slab_truncation', 'spherical_truncation'],
                    ['use_xdat', 'dont_use_vxcdat'], 
                    ['ecdel', 'cvfit'], ['ecdel_outer', 'cvfit_outer'],
                    ['use_symmetries_q_grid', 'no_symmetries_q_grid'],
                    ['die_outside_sphere', 'ignore_outside_sphere']]
                },
            'kernel': {
                'int': ['number_val_bands', 'number_cond_bands'],
                'float': ['screened_coulomb_cutoff', 'bare_coulomb_cutoff',
                    'fermi_level', 'coulomb_truncation_radius'],
                'str': ['partial_blocks'],
                'bool': ['fermi_level_absolute', 'fermi_level_relative', 
                    'screening_metal', 'screening_graphene',
                    'screening_semiconductor', 'cell_box_truncation',
                    'cell_wire_truncation', 'cell_slab_truncation',
                    'spherical_truncation', 'use_symmetries_coarse_grid',
                    'no_symmetries_coarse_grid', 'write_vcoul', 'comm_mpi',
                    'comm_disk', 'low_comm', 'low_memory', 'fullbz_replace', 
                    'fullbz_write', 'die_outside_sphere', 
                    'ignore_outside_sphere', 'read_kpoints'],
                'exclusive': [['fermi_level_absolute', 'fermi_level_relative'],
                    ['screening_metal', 'screening_graphene', 
                    'screening_semiconductor'], ['cell_box_truncation',
                    'cell_wire_truncation', 'cell_slab_truncation',
                    'spherical_truncation'], ['use_symmetries_coarse_grid',
                    'no_symmetries_coarse_grid'], ['comm_mpi', 'comm_disk'],
                    ['die_outside_sphere', 'ignore_outside_sphere']]},
            'absorption': {
                'int': ['number_val_bands_fine', 'number_val_bands_coarse',
                    'number_cond_bands_fine', 'number_cond_bands_coarse',
                    'lowest_occupied_band', 'highest_occupied_band',
                    'coarse_grid_points', 'number_eigenvalues', 
                    'number_iterations', 'write_eigenvectors'],
                'float': ['fermi_level', 'evs', 'ev0', 'evdel', 'ecs', 'ec0',
                    'ecdel', 'cvfit', 'no_symmetries_fine_grid', 
                    'use_symmetries_fine_grid', 'no_symmetries_shifted_grid',
                    'use_symmetries_shifted_grid', 'no_symmetries_coarse_grid',
                    'use_symmetries_coarse_grid', 'coulomb_truncation_radius',
                    'cell_average_cutoff', 'energy_resolution', 
                    'energy_resolution_sigma', 'energy_resolution_gamma', 
                    'kernel_scaling'],
                'str': ['regular_grid', 'q_shift', 'polarization'],
                'bool': ['fermi_level_absolute', 'fermi_level_relative',
                    'no_symmetries_fine_grid', 'use_symmetries_fine_grid',
                    'no_symmetries_shifted_grid', 'use_symmetries_shifted_grid',
                    'no_symmetries_coarse_grid', 'use_symmetries_coarse_grid',
                    'comm_mpi', 'comm_disk', 'eqp_corrections', 
                    'eqp_co_corrections', 'fullbz_replace', 'fullbz_write', 
                    'read_kpoints', 'diagonalization', 'haydock',
                    'screening_metal', 'screening_graphene', 
                    'screening_semiconductor', 'cell_box_truncation', 
                    'cell_wire_truncation', 'cell_slab_truncation',
                    'spherical_truncation', 'use_velocity', 'use_momentum', 
                    'use_dos', 'read_vmtxel', 'read_eps2_moments', 
                    'read_eigenvalues', 'noeh_only', 'read_dtmat', 
                    'voigt_broadening', 'read_epsdiag', 'spin_triplet', 
                    'spin_singlet', 'local_fields', 'write_vcoul',
                    'skip_interpolation', 'degeneracy_check_override', 
                    'average_w', 'lorentzian_broadening', 
                    'gaussian_broadening'], 
                'exclusive': [['fermi_level_absolute', 'fermi_level_relative'],
                    ['ecdel', 'cvfit'], ['no_symmetries_fine_grid',
                    'use_symmetries_fine_grid'], ['no_symmetries_shifted_grid',
                    'use_symmetries_shifted_grid'], 
                    ['no_symmetries_coarse_grid', 'use_symmetries_coarse_grid'],
                    ['diagonalization', 'haydock'], ['screening_metal',
                    'screening_graphene', 'screening_semiconductor'], 
                    ['cell_box_truncation', 'cell_wire_truncation', 
                    'cell_slab_truncation', 'spherical_truncation'],
                    ['use_velocity', 'use_momentum', 'use_dos'],
                    ['spin_triplet', 'spin_singlet', 'local_fields'],
                    ['lorentzian_broadening', 'gaussian_broadening']]
                }}

        run_dict = valid_params[self.run_type]
        key_found = False
        
        def check_exclusives_and_set(k,v):
            rd2 = cp.deepcopy(run_dict['exclusive'])
            # Check to see if Key is mutually exclusive with other keys
            for i,l in enumerate(rd2):
                if k in l:
                    l.remove(k)

                    # With Key removed from list, check if other values
                    # are already used inside Input Set Parameters
                    for j in l:
                        if isp_dict.has_key(j):
                            raise Exception('Parameters "{}" & "{}" '.format(
                                    k, j) + 'cannot be used together.  ' +
                                    'They are mutually exclusive. ') 

            # Set Key/Value pair if conflicts were not found
            isp_dict.update({k:v})

        for i in run_dict.keys():
            if key in run_dict[i]:
                key_found = True
                key_dict = i
        if not key_found:
            ignore_list = ['kpoints', 'qshift', 
                    'isp', 'run_type', 'mat_type', 'filename',  
                    'kps', 'kgrid', 'cmplx_real','qemf_dir',
                    'config_file','qe_dirs', 'bgw_dirs', 'reduce_structure']
            if key in ignore_list:
                self.__dict__['params'].update({key:val})
            if key not in ignore_list:
                print("\nKey: '{}' not found in valid parameters. ".format(key)+
                        "Setting Key/Val pair as class attribute.")
            self.__dict__.update({key: val})

        else:
            if key_dict == 'int':
                try:
                    check_exclusives_and_set(key, int(val))
                except:
                    raise TypeError("Value for Key: {} ".format(key)+
                            "should be of Type 'int', not {}.".format(
                                type(val)))
            elif key_dict == 'float':
                try:
                    print "gk: calling check_exclusives_and_set with key=",key,"val=",val
                    check_exclusives_and_set(key, float(val))
                except:
                    raise TypeError("Value for Key: {} ".format(key)+
                            "should be of Type 'float', not {}.".format(
                                type(val)))
            elif key_dict == 'bool':
                if not isinstance(val,bool):
                    raise TypeError("Value for Key: {} ".format(key)+
                            "should be of Type 'bool', not {}.".format(
                                type(val)))
                check_exclusives_and_set(key, '')
            else:
                check_exclusives_and_set(key, val)
        
    @staticmethod
    def from_file(file):
        print "gk: in BgwInput.from_file, file=",file
        with zopen(file, 'rt') as f:
            return BgwInput.from_string(f.read(), file)

    @staticmethod
    def from_string(string, filename):
        """
        Reads BgwInput Object from string.

        Args:
            string (str):   BerkeleyGW Input String
            filename (str): Output filename

        Returns:
            BgwInput object
        """
        print "gk: in BgwInput.from_string"
        run_type = filename.split('/')[-1].split('.')[0]
        lines = list(string.splitlines())
        kps = None
        isp_params = {}
        for line in lines:
            m = re.match("(\w+)\s*(.*)", line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).strip() if m.group(2) else ''

                if key != 'begin' and key != 'end':
                    isp_params[key] = val
                elif key == 'begin':
                    kps = lines.index(line) + 1
                elif key == 'end':
                    kp_end = lines.index(line)

        if kps:
            kps = lines[kps:kp_end]

        params = {'isp': isp_params, 
                'kps': kps,  
                'filename': "{}.inp".format(run_type)}

        return BgwInput(**params)

    @staticmethod
    def from_directory(d):
        in_files = glob.glob("{}/*.inp".format(d))
        if len(in_files) > 1:
            print("Found more than one Input file, using {}".format(
                    in_files[0]))
        return BgwInput.from_file(in_files[0])


