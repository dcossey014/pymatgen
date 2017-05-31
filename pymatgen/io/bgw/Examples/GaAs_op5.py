#!/usr/bin/env python

from fireworks import Firework, Workflow

import os, copy
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.bgw.interfaces import QeMeanFieldTask, BgwAbsTask, BgwFirework, BgwWorkflow
from pymatgen.io.bgw.inputs import BgwInput
from pymatgen.io.bgw.pwscf_tasks import BgwCustodianTask

s = Structure.from_file('GaAs_sg216_exp.cif')
finder = SpacegroupAnalyzer(s)

s_prim = finder.get_primitive_standard_structure()
#s_prim.translate_sites([0,1], [-0.125,-0.125,-0.125], to_unit_cell=False)
print s_prim
print s_prim.__repr__()

# get rid of spaces for pwscf and pw2bgw file prefixes
formula_prefix=s_prim.formula.replace(' ','')

kpoints = [6,6,6]
qshift=[0.0, 0.0, 0.001]
fftw_grid = [24,24,24]

pseudo_dir = '/app/ccm/PPs/ONCVPSP/QE'
mpi_cmd = 'mpiexec_mpt -n 36'

pw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw.x'
pw2bgw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw2bgw.x'

cmplx_real='cmplx'

mean_field_tasks={'qe_tasks' : ['scf','wfn','wfnq','wfn_co','wfn_fi','wfnq_fi'] }

scf_control={"calculation" : "scf", 
                "pseudo_dir": pseudo_dir,
                'prefix': formula_prefix,
                'tstress': True, 
                'tprnfor': True,
                'verbosity': 'high'}

wfn_control={"calculation" : "bands",
                "pseudo_dir": pseudo_dir,
                'prefix': formula_prefix,
                'tstress': True, 
                'tprnfor': True}

qemft_control={}

qemft_control={}
qemft_control['scf']=scf_control

for i in ['wfn', 'wfnq', 'wfn_co', 'wfn_fi', 'wfnq_fi']:
    qemft_control[i] = copy.deepcopy(wfn_control)

scf_system={'ecutwfc': 60.0 }

qemft_system = {}
for i in ['scf', 'wfn', 'wfnq', 'wfn_co', 'wfn_fi', 'wfnq_fi']:
    qemft_system[i] = copy.deepcopy(scf_system)

qemft_system['wfn']['nbnd'] = 40
qemft_system['wfnq']['nbnd'] = 16
qemft_system['wfn_co']['nbnd'] = 45
qemft_system['wfn_fi']['nbnd'] = 16
qemft_system['wfnq_fi']['nbnd'] = 16


scf_electrons={'conv_thr': 1.0e-10, 'diago_full_acc': True, 'diagonalization': 'cg'}
qemft_electrons={}
for i in ['scf', 'wfn', 'wfnq', 'wfn_co', 'wfn_fi', 'wfnq_fi']:
    qemft_electrons[i] = copy.deepcopy(scf_electrons)

wfn_file='wfn.'+cmplx_real
if cmplx_real == 'real':
    real_or_complex=1
elif cmplx_real == 'cmplx':
    real_or_complex=2
else:
    print "unknown cmplx_real"
    exit()

basic_in_pw2bgw = {
    'prefix':  formula_prefix,
    'wfng_flag': True, 
    'real_or_complex': real_or_complex,  
    'wfng_file': wfn_file,
    'wfng_kgrid': False
}

# Set PostProcessing QE Dictionary Parameters
qemft_pw2bgw={}
for i in ['wfn', 'wfnq', 'wfn_co', 'wfn_fi', 'wfnq_fi']:
    qemft_pw2bgw[i] = copy.deepcopy(basic_in_pw2bgw)

rhog_file='rho.'+cmplx_real

qe_wfn_co_in_pw2bgw = qemft_pw2bgw['wfn_co']
qe_wfn_co_in_pw2bgw['rhog_flag'] = True
qe_wfn_co_in_pw2bgw['rhog_file']=rhog_file
qe_wfn_co_in_pw2bgw['vxc_flag']=True
qe_wfn_co_in_pw2bgw['vxc_file']='vxc.dat'
qe_wfn_co_in_pw2bgw['vxc_diag_nmin']=1
qe_wfn_co_in_pw2bgw['vxc_diag_nmax']=45
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmin']=0
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmax']=0

qemft = QeMeanFieldTask(structure=s_prim, kpoints=kpoints, pseudo_dir=pseudo_dir,
        mpi_cmd=mpi_cmd, pw_cmd=pw_cmd, pw2bgw_cmd=pw2bgw_cmd, qshift=qshift,
        fftw_grid=fftw_grid, qe_control=qemft_control, qe_system=qemft_system,
        qe_electrons=qemft_electrons, qe_pw2bgw=qemft_pw2bgw, mf_tasks=mean_field_tasks)

qemf_fw = Firework([qemft], name="QeMeanFieldExample")

#qemft.write_inputs()
bgw_wf = BgwWorkflow(qemf_fw, name="QETask")

#bgw_wf.to_file('qe_gaas.yaml')

########################
#Begin BGW FireWorks   #
########################

#######################
# Epsilon Calculation #
#######################

if cmplx_real == 'cmplx':
    cmplx_bool=True
else:
    cmplx_bool=False

# set to None if you don't want to use this
prev_qemf_dir=None
#prev_qemf_dir='/home/dec014/OP4/ESPRESSO'

eps_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
                kpoints=kpoints, qshift=qshift, filename='epsilon.inp',
                qemf_dir=prev_qemf_dir)

eps_inp.epsilon_cutoff=20.0
#eps_inp.number_bands = 39

eps_fw = BgwFirework(eps_inp, name="Epsilon Task", complex=cmplx_bool)

sig_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
                kpoints=kpoints, qshift=qshift, filename='sigma.inp')

sig_inp.screened_coulomb_cutoff = 10.0
sig_inp.bare_coulomb_cutoff = 25.0
sig_inp.band_index_min = 1
sig_inp.screening_semiconductor = True

# what is this?  It looks out of place 
bgw_dir = "/apps/ccm/opt/BerkeleyGW-1.1-beta2/bin/"
eqp = bgw_dir+'eqp.py'
ppx = ' '.join([eqp, 'eqp1', './sigma_hp.log', './eqp_co.dat'])

sig_fw = BgwFirework(sig_inp, name="Sigma Task", ppx=ppx, complex=cmplx_bool)

krn_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, filename='kernel.inp', cmplx_real=cmplx_real)

krn_inp.screened_coulomb_cutoff = 10.0
krn_inp.bare_coulomb_cutoff = 25.0
krn_inp.use_symmetries_coarse_grid = True
krn_inp.screening_semiconductor = True
krn_inp.number_cond_bands = 10
krn_inp.number_val_bands = 4

krn_fw = BgwFirework(krn_inp, name="Kernel Task", complex=cmplx_bool, mpi_cmd='mpiexec_mpt -n 72')
krn_fw.add_spec('_queueadapter', {'nnodes': 2})

abs_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, filename='absorption.inp', cmplx_real=cmplx_real)

abs_inp.diagonalization = True
abs_inp.use_symmetries_coarse_grid = True
abs_inp.no_symmetries_fine_grid = True
abs_inp.no_symmetries_shifted_grid = True
abs_inp.screening_semiconductor = True
abs_inp.use_velocity = True
abs_inp.gaussian_broadening = True
abs_inp.energy_resolution = 0.15
abs_inp.eqp_co_corrections = True
abs_inp.number_val_bands_fine = 4
abs_inp.number_val_bands_coarse = 4
abs_inp.number_cond_bands_fine = 6
abs_inp.number_cond_bands_coarse = 10

abs_fw = BgwFirework(abs_inp, name='Absorption Task', complex=True, mpi_cmd='mpiexec_mpt -n 72')
abs_fw.add_spec('_queueadapter', {'walltime': '4:00:00', 'nnodes': 4, 'ppnode': 18,
                                'mppwidth': 18, 'queue': 'standard'})

bgw_wf2 = BgwWorkflow(qemf_fw, eps_fw, sig_fw, krn_fw, abs_fw, 
                    name="QE/BGW Optical Task")

#bgw_wf2 = BgwWorkflow(eps_fw, sig_fw, krn_fw, abs_fw,
#                    name="QE/BGW Optical Task")

#bgw_wf2 = BgwWorkflow(qemf_fw, eps_fw,
#                    name="QE/BGW Optical Task")

bgw_wf2.to_file('gaas_op5.yaml')
