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

from pymatgen.io.bgw.kgrid import generate_kpath

from monty.io import zopen
from monty.dev import deprecated
from monty.json import MSONable
from monty.serialization import loadfn
from collections import defaultdict, OrderedDict
from pymatgen.io.bgw import Kgrid
from fireworks import FireTaskBase, explicit_serialize, FWAction
from pymatgen import Structure


@explicit_serialize
class BgwInputTask(FireTaskBase):
    '''
    Fireworks Task for writing BerkeleyGW Input Files.

    Required Parameters:
        structure           :   An input structure in pymatgen's Structure format
        pseudo_dir  (str)   :   Directory where PseudoPotential Files are stored
        input_set_params (dict):Dictionary of input set parameters for the BerkeleyGW
                                job to be run
        out_file    (str)   :   Filename to write Input file to

    Optional Parameters:
        kpoints (array)     :   List of Kpoints used in Kgrid.x to create WFN.
                                This is required for Epsilon and Sigma calculations.
        qshift  (array)     :   Qshift used in Kgrid.x for WFN Qshifted grids.
                                This is required for Epsilon calculations.
        mat_type (str)      :   Options: 'metal' or 'semiconductor'.  This is required
                                for Epsilon calculations.

    '''
    
    required_params = ['structure', 'pseudo_dir', 'input_set_params', 'out_file']
    optional_params = ['kpoints', 'qshift', 'mat_type', 'cmplx_real', 'qemf_dir', 'config_file']

    def __init__(self, params): 
        self.structure = Structure.from_dict(params.get('structure').as_dict()) if params.get('structure', '') else ''
        self.pseudo_dir = params.get('pseudo_dir')
        self.kpoints = params.get('kpoints', None)
        self.qshift = params.get('qshift', None)
        self.isp = params.get('input_set_params')
        self.run_type = params.get('run_type', None)
        self.mat_type = params.get('mat_type', 'metal')
        self.filename = params.get('out_file')
        self.cmplx_real = params.get('cmplx_real')
        self.qemf_dir = params.get('qemf_dir')
        self.kps = params.get('kps', None)
        #gk: why is occupied_bands set here?
        self.occupied_bands = params.get('occupied_bands', 0)
        self.config_file = params.get('config_file')

        #print "gk: in BgwInputTask.__init__ from params, self.occupied_bands = ",self.occupied_bands

        #Build Pseudo Dictionary and List of PPs files if given Structure
        if self.structure:
            self.pseudo = {}
            self.pseudo_files = []
            for i in self.structure.symbol_set: 
                pattern = "{}*.upf".format(i)
                pps = []
                for j in os.listdir(self.pseudo_dir):
                    m = re.match(re.compile(fnmatch.translate(pattern), 
                                    re.IGNORECASE), j)
                    if m:
                        pps.append(os.path.join(self.pseudo_dir,m.group(0)))
                self.pseudo_files.append(pps[-1])
                self.pseudo[i.encode('ascii', 'ignore')] = pps[-1].split('/')[-1]

        params = {'structure': self.structure, 'pseudo_dir': self.pseudo_dir,
                'kpoints': self.kpoints, 'qshift': self.qshift, 
                'input_set_params': self.isp, 'run_type': self.run_type,
                'mat_type':self.mat_type, 'out_file': self.filename,
                'occupied_bands': self.occupied_bands, 
                'kps': self.kps, 'cmplx_real': self.cmplx_real,'qemf_dir': self.qemf_dir}
        #super(BgwInputTask, self).__init__()
        #self.check_write_params()
        self.update(params)
        print "gk: in BgwInput.__init__(), self.occupied_bands=",self.occupied_bands
        
    def __str__(self):
        out = []
        def write_kpoints():
            if not self.kps:

                '''
                #print("Kpoints: %s" %(self.kpoints))
                self.kgrid = Kgrid(self.structure, kpoints=self.kpoints, 
                                offset_type="gamma")
                self.kps = self.kgrid.generate_kpoints('WFN')
                drop_line = self.kps.pop(0)
    
                qtype = 2 if type == 'metal' else 1
            
                for i,j in enumerate(self.kps):
                    k = j.split()
                    self.kps[i] = "{0:.6f}  {1:.6f}  {2:.6f}  1.0\n".format(float(k[0]),
                                        float(k[1]), float(k[2]))
            
                if self.run_type == 'epsilon':
                    self.kps[0] = "{0:.6f}  {1:.6f}  {2:.6f}  1.0  {3}\n".format(
                                    self.qshift[0], self.qshift[1], self.qshift[2],
                                    qtype)

                    for i,j in enumerate(self.kps[1:], start=1):
                        k = j.split()
                        self.kps[i] = "{0:.6f}  {1:.6f}  {2:.6f}  {3:>3.1f}  0\n".format(
                                        float(k[0]), float(k[1]), float(k[2]), float(k[3]))
            '''
            for kp in self.kps:
                out.append("  {}".format(kp.strip()))
            out.append('end')

        def write_band_occ():
            print "gk in BgwInput.__str__.write_band_occ(), self.occupied_bands=",self.occupied_bands
            out.append("band_occupation {}*1 {}*0".format(
                    int(self.occupied_bands), 
                    int(self.isp['number_bands']) - int(self.occupied_bands)))

        for k1 in self.isp.keys():
            out.append("{}  {}".format(k1, self.isp[k1]))
        print "gk: first out = ", out

        if 'epsilon' in self.run_type:
            if self.occupied_bands != None:
                write_band_occ()
            out.append('begin qpoints')
            write_kpoints()
        elif 'sigma' in self.run_type:
            if self.occupied_bands != None:
                write_band_occ()
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
        self.check_params()
        with open(filename, 'w') as f:
            f.write(self.__str__()+'\n')

    def check_params(self):
        #gk: probably none of these are required in BGW v 1.2
        eps_params = ['epsilon_cutoff', 'number_bands']
        #gk: probably none of these are required
        #sig_params = ['screened_coulomb_cutoff', 'bare_coulomb_cutoff',
        #            'number_bands', 'band_index_min', 'band_index_max']
        sig_params = ['number_bands', 'band_index_min', 'band_index_max']
        krn_params = ['number_val_bands', 'number_cond_bands', 
                    'screened_coulomb_cutoff', 'bare_coulomb_cutoff']
        abs_params = ['number_val_bands_coarse', 'number_val_bands_fine',
                    'number_cond_bands_coarse', 'number_cond_bands_fine',
                    'energy_resolution']
        sig2wan_params = ['sigma_out', 'wannier_input', 'spin_component',
                    'eqp_num', 'num_bands']

        if 'epsilon' in self.run_type:
            params = list(eps_params)
        elif 'sigma' in self.run_type:
            #gk: sigma picks up number_bands from epsilon? 
            params = list(sig_params)
            self.isp['band_index_min'] = self.isp.get('band_index_min',1)
            #self.isp['band_index_max'] = self.isp.get('band_index_max', 
            #                            self.isp.get('number_bands', 0))
            bim=self.isp.get('band_index_max',self.isp.get('number_bands', 0))
            bm=self.isp.get('number_bands')
            bim=bm if bim>bm else bim
            self.isp['band_index_max']=bim
            print "gk: in BgwInput.check params"
            print "gk:     band_index_min=", self.isp.get('band_index_min')
            print "gk:     band_index_max=", self.isp.get('band_index_max')
            print "gk:     number_bands=", self.isp.get('number_bands')
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

        #with open(filename, 'w') as f:
        #    f.write(self.__str__()+'\n')

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
            #wfn = self.qe_dirs['wfn']
            force_link(os.path.join(wfn_co, 'vxc.dat'),
                    './vxc.dat')
            rho_file='rho.'+self.cmplx_real
            force_link(os.path.join(wfn_co, rho_file),
                    './RHO')
            #gk: number_bands from here?
            #force_link(os.path.join(wfn_co, wfn_co_file),
            #        './WFN_outer')
            #gk: or number_bands from here?
            #gk: WFN_inner could be wfn_co or wfn....????
            force_link(os.path.join(wfn_co, wfn_file), 
                    './WFN_inner')
            extra_links('eps')
            #gk: number_bands checked and adjusted here
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

            # Old Code
            '''
            matched_band_index=0
            matched_band_diff=9999
            matched_band_sign=0
            for band_index in allowed:
                diff=user_spec-band_index
                sign=1
                if diff < 0: sign=-1
                if abs(diff) < matched_band_diff:      
                    matched_band_diff=abs(diff)
                    matched_band_sign=sign
                    matched_band_index=band_index
                elif abs(diff) == matched_band_diff:
                    if sign > matched_band_sign:
                        matched_band_diff=abs(diff)
                        matched_band_sign=sign
                        matched_band_index=band_index
            return matched_band_index
            '''

        if not self.occupied_bands:
            self.occupied_bands = allowed_vbands[-1]
            print "occupied_bands not set.  Setting including all occupied bands, occupied_bands=",self.occupied_bands
        else:
            print "User specified occupied_bands =", self.occupied_bands
            print "Checking for allowed bands closest to occupied_bands."
            matched_vband = match_bands(self.occupied_bands,allowed_vbands)
            if matched_vband == self.occupied_bands:
                print "User specified occupied_bands is degeneracy allowed"
            else:
                print "occupied_bands is reset to to degeneracy allowed value,",matched_vband
                self.occupied_bands=matched_vband
    
        if 'epsilon' in self.filename:
            user_nbands=int(self.isp.get('number_bands', 0))
            print "For Epsilon user specified number_bands =", user_nbands
            if user_nbands < self.occupied_bands:
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
            if user_nbands < self.occupied_bands:
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
            if max_nbands < self.occupied_bands:
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
        
        
    def dry_run(self):

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
        structure (Structure):  Pymatgen Structure to be used in the BerkeleyGW
                                calculation.
        pseudo_dir (str):       Directory where PseudoPotential files are
                                stored.
        filename (str):         Filename to write Input file to.
        kpoints (array):        List of Kpoints used in Kgrid.x to create WFN.
                                This is required for Epsilon and Sigma calculations.
        qshift  (array):        Qshift used in Kgrid.x for WFN Qshifted grids.
                                This is required for Epsilon calculations.
        mat_type (str):         Options: 'metal' or 'semiconductor'.  This is
                                required for Epsilon calculations.
                                Default: semiconductor
    '''

    def __init__(self, structure, pseudo_dir, isp={}, cmplx_real='cmplx',
                kpoints=None, qshift=None, mat_type='semiconductor',
                kps=None, occupied_bands = None, filename=None,
                qemf_dir=None):

        self.__dict__['isp'] = isp if isp else {}
        self.__dict__['run_type'] = filename.split('/')[-1].split('.')[0]
        self.__dict__['params'] = {'structure': structure, 'pseudo_dir': pseudo_dir,
                'kpoints': kpoints, 'qshift': qshift, 'kps': kps,
                'input_set_params': self.isp, 'mat_type': mat_type,
                'out_file': filename, 'run_type': self.run_type,
                'occupied_bands': occupied_bands, 'cmplx_real': cmplx_real,
                'qemf_dir': qemf_dir}

        super(BgwInput, self).__init__(self.params)



    def write_control_file(self):
        filename = self.params.get('out_file')
        self.write_file(filename)

    def __setitem__(self, key, val):
        self.proc_key_val(key.strip(), val.strip()) if isinstance(
                val, six.string_types) else self.proc_key_val(key.strip(), val)

    def __setattr__(self, key, val):
        self.proc_key_val(key.strip(), val.strip()) if isinstance(
                val, six.string_types) else self.proc_key_val(key.strip(), val)

    def proc_key_val(self, key, val):
        #print("key: {}\nVal: {}\n\n".format(key, val))
        isp_dict = self.params.get('input_set_params', {})
        #print("Input Dictionary:  {}\n\n".format(isp_dict))

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
            ignore_list = ['structure', 'pseudo_dir', 'kpoints', 'qshift', 
                    'isp', 'run_type', 'mat_type', 'filename', 'pseudo', 
                    'pseudo_files', 'occupied_bands', 'kps', 'kgrid', 'cmplx_real','qemf_dir',
                    'config_file','qe_dirs', 'bgw_dirs']
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
                    check_exclusives_and_set(key, float(val))
                except:
                    raise TypeError("Value for Key: {} ".format(key)+
                            "should be of Type 'float', not {}.".format(
                                type(val)))
            elif key_dict == 'bool':
                if type(val) != bool:
                    raise TypeError("Value for Key: {} ".format(key)+
                            "should be of Type 'bool', not {}.".format(
                                type(val)))
                check_exclusives_and_set(key, '')
            else:
                check_exclusives_and_set(key, val)
            
        
    @staticmethod
    def from_file(file):
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
        run_type = filename.split('/')[-1].split('.')[0]
        lines = list(string.splitlines())
        kps = None
        occupied_bands = None
        isp_params = {}
        for line in lines:
            m = re.match("(\w+)\s*(.*)", line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).strip() if m.group(2) else ''

                if 'band_occupation' in key:
                    n = re.match("(\w+)\*.*", val)
                    occupied_bands = n.group(1).strip()
                elif key != 'begin' and key != 'end':
                    isp_params[key] = val
                elif key == 'begin':
                    kps = lines.index(line) + 1
                elif key == 'end':
                    kp_end = lines.index(line)

        if kps:
            kps = lines[kps:kp_end]

        params = {'pseudo_dir': '', 'isp': isp_params, 
                'kps': kps, 'occupied_bands': occupied_bands, 
                'filename': "{}.inp".format(run_type)}

        return BgwInput('', **params)

    @staticmethod
    def from_directory(d):
        in_files = glob.glob("{}/*.inp".format(d))
        if len(in_files) > 1:
            print("Found more than one Input file, using {}".format(
                    in_files[0]))
        return BgwInput.from_file(in_files[0])


@deprecated(BgwInput, 'Use BgwInput instead')
class BGWInput():
    def __init__(self, structure, pseudo_dir, input_set_params={},
                run_type=None, kpoints=None, qshift=None, type='metal',
                filename=None):

        BgwInput.__init__(self, structure, pseudo_dir, input_set_params={},
                run_type=None, kpoints=None, qshift=None, type='metal',
                filename=out_file)


class PwPpInput(abc.ABCMeta, MSONable):
    #TODO   Look into MRO and inheritance vs python version
    #       maybe use six for portability across versions like VASP?
    '''
    Docstring
    '''

    def as_dict(self):
        d = MSONable.as_dict(self)
        return d

    @abc.abstractmethod
    def write_input(self, filename):
        """Write Input Files"""
        pass

class PwBandsInput(PwPpInput):
    """
    Docstring
    """
    def __init__(self, bands=None):
        self.bands = bands


    def __str__(self):
        return

    def write_input(self, filename):
        pass

class PW2BGWInput(object):
    '''
    Base Input file Class for PW2BGW post-processing.
    '''

    def __init__(self, structure, pw2bgw_input=None, kpoints=None, k_offset=None, qshift=None ):
        '''
        Initializes a PW2BGW input file.

        Args:
            structure (Structure): Input Stucture.
            pw2bgw_input (dict): Input parameters.  Refer to official PW2BGW 
                Quantum Espresso doc for supported parameters. Defaults to
                {'prefix': Structure.formula, 'wfng_flag': True,
                'wfng_kgrid': False}
        '''

        self.structure = structure
        self.sections = {}
        self.sections['input_pw2bgw'] = pw2bgw_input or {'prefix': Structure.formula,
                                'wfng_flag': True, 'wfng_kgrid': False}
        self.kpoints=kpoints

        if self.kpoints:
            self.sections['input_pw2bgw']['wfng_kgrid'] = True
            self.sections['input_pw2bgw']['wfng_nk1'] = int(self.kpoints[0])
            self.sections['input_pw2bgw']['wfng_nk2'] = int(self.kpoints[1])
            self.sections['input_pw2bgw']['wfng_nk3'] = int(self.kpoints[2])

            self.k_offset = k_offset if k_offset else "Monkhorst-Pack"
            self.qshift = qshift if qshift else [0, 0, 0]

            if isinstance(self.k_offset, list):
                self.sections['input_pw2bgw']['wfng_dk1'] = self.k_offset[0]
                self.sections['input_pw2bgw']['wfng_dk2'] = self.k_offset[1]
                self.sections['input_pw2bgw']['wfng_dk3'] = self.k_offset[2]
            elif "monkhorst" in self.k_offset.lower(): 
                self.sections['input_pw2bgw']['wfng_dk1'] = 0.5 + self.qshift[0]*self.kpoints[0]
                self.sections['input_pw2bgw']['wfng_dk2'] = 0.5 + self.qshift[1]*self.kpoints[1]
                self.sections['input_pw2bgw']['wfng_dk3'] = 0.5 + self.qshift[2]*self.kpoints[2]
            elif 'random' in self.k_offset.lower():
                self.sections['input_pw2bgw']['wfng_dk1'] = 0.47 + self.qshift[0]*self.kpoints[0]
                self.sections['input_pw2bgw']['wfng_dk2'] = 0.37 + self.qshift[1]*self.kpoints[1]
                self.sections['input_pw2bgw']['wfng_dk3'] = 0.32 + self.qshift[2]*self.kpoints[2]
            else:
                self.sections['input_pw2bgw']['wfng_dk1'] = self.qshift[0]*self.kpoints[0]
                self.sections['input_pw2bgw']['wfng_dk2'] = self.qshift[1]*self.kpoints[1]
                self.sections['input_pw2bgw']['wfng_dk3'] = self.qshift[2]*self.kpoints[2]

    def __str__(self):
        out = []
        def to_str(v):
            if isinstance(v, six.string_types):
                return "'%s'" % v
            return v

        for k1 in ['input_pw2bgw']:
            v1 = self.sections[k1]
            out.append("&%s" % k1.upper())
            sub = []
            for k2 in sorted(v1.keys()):
                sub.append("   %s = %s" % (k2, to_str(v1[k2])))
            sub.append("/")
            out.append(",\n".join(sub))
        return "\n".join(out)

    def write_input(self, filename):
        '''
        Write the PW2BGW input file.

        Args:
            filename (str): The string filename to output to.
        '''
        with open(filename, 'w') as fout:
            fout.write(self.__str__()+"\n")


