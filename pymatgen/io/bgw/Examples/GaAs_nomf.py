#!/usr/bin/env python

from fireworks import Firework, Workflow

import os, copy
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.bgw.interfaces import BgwFirework, BgwWorkflow
from pymatgen.io.espresso.interfaces import QeMeanFieldTask
from pymatgen.io.bgw.inputs import BgwInput
from pymatgen.io.bgw.custodian_jobs import BgwCustodianTask

s = Structure.from_file('GaAs_sg216_exp.cif')
finder = SpacegroupAnalyzer(s)

s_prim = finder.get_primitive_standard_structure()
#s_prim.translate_sites([0,1], [-0.125,-0.125,-0.125], to_unit_cell=False)
print s_prim
print s_prim.__repr__()

# get rid of spaces for pwscf and pw2bgw file prefixes
formula_prefix=s_prim.formula.replace(' ','')

kpoints = [5,5,5]
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
                'tprnfor': True}

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

scf_system={'ecutwfc': 25.0 }

qemft_system = {}
for i in ['scf', 'wfn', 'wfnq', 'wfn_co', 'wfn_fi', 'wfnq_fi']:
    qemft_system[i] = copy.deepcopy(scf_system)

qemft_system['wfn']['nbnd'] = 30
qemft_system['wfnq']['nbnd'] = 7
qemft_system['wfn_co']['nbnd'] = 33
qemft_system['wfn_fi']['nbnd'] = 11
qemft_system['wfnq_fi']['nbnd'] = 4


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
qe_wfn_co_in_pw2bgw['vxc_diag_nmax']=14
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmin']=0
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmax']=0

#qemft = QeMeanFieldTask(structure=s_prim, kpoints=kpoints, pseudo_dir=pseudo_dir,
#        mpi_cmd=mpi_cmd, pw_cmd=pw_cmd, pw2bgw_cmd=pw2bgw_cmd, qshift=qshift,
#        fftw_grid=fftw_grid, qe_control=qemft_control, qe_system=qemft_system,
#        qe_electrons=qemft_electrons, qe_pw2bgw=qemft_pw2bgw, mf_tasks=mean_field_tasks)
#
#qemf_fw = Firework([qemft], name="QeMeanFieldExample")

#qemft.write_inputs()
#bgw_wf = BgwWorkflow(qemf_fw, name="QETask")
#
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
prev_qemf_dir='/workspace/gkedz/BGWFire/GaAs/QE4/ESPRESSO'

eps_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
                kpoints=kpoints, qshift=qshift, filename='epsilon.inp',
                qemf_dir=prev_qemf_dir)

eps_inp.epsilon_cutoff=20.0
eps_inp.number_bands = 27

# use this to test
#eps_inp.dry_run()


eps_fw = BgwFirework(eps_inp, name="Epsilon Task", complex=cmplx_bool, mpi_cmd='mpiexec_mpt -n 36')

bgw_wf2 = BgwWorkflow(eps_fw, name="QE/BGW Optical Task")
bgw_wf2.to_file('gaas_nomf.yaml')