#!/usr/bin/env python

from fireworks import Firework, Workflow

import os, copy, string, json, pprint
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.bgw.interfaces import QeMeanFieldTask, BgwAbsTask, BgwFirework, BgwWorkflow
from pymatgen.io.bgw.inputs import BgwInput

s = Structure.from_file('GaAs_sg216_exp.cif')
finder = SpacegroupAnalyzer(s)

s_prim = finder.get_primitive_standard_structure()
#s_prim.translate_sites([0,1], [-0.125,-0.125,-0.125], to_unit_cell=False)
print s_prim
print s_prim.__repr__()

# get rid of spaces for pwscf and pw2bgw file prefixes
formula_prefix=s_prim.formula.replace(' ','')

kpts = [12,12,12]

kpoints = {}
kpoints['scf'] = kpts
kpoints['wfn'] = [24,24,24]

qshift=[0.0, 0.0, 0.001]
fftw_grid = [24,24,24]

pseudo_dir = '/app/ccm/PPs/ONCVPSP/QE'
mpi_cmd = 'mpiexec_mpt -n 36'

pw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw.x'
pw2bgw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw2bgw.x'

cmplx_real='cmplx'

mean_field_tasks={'qe_tasks' : ['scf','wfn'] }

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

for i in ['wfn']:
    qemft_control[i] = copy.deepcopy(wfn_control)

scf_system={'ecutwfc': 60.0 }

qemft_system = {}
for i in ['scf', 'wfn']:
    qemft_system[i] = copy.deepcopy(scf_system)

qemft_system['wfn']['nbnd'] = 40

scf_electrons={'conv_thr': 1.0e-10, 'diago_full_acc': True, 'diagonalization': 'cg'}
qemft_electrons={}
for i in ['scf', 'wfn']:
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
for i in ['wfn']:
    qemft_pw2bgw[i] = copy.deepcopy(basic_in_pw2bgw)

qemft = QeMeanFieldTask(structure=s_prim, kpoints=kpoints, pseudo_dir=pseudo_dir,
        mpi_cmd=mpi_cmd, pw_cmd=pw_cmd, pw2bgw_cmd=pw2bgw_cmd, qshift=qshift,
        fftw_grid=fftw_grid, qe_control=qemft_control, qe_system=qemft_system,
        qe_electrons=qemft_electrons, qe_pw2bgw=qemft_pw2bgw, 
        krgrid_offset_type='Gamma', mf_tasks=mean_field_tasks,
        bandstructure_kpoint_path=True, num_kpoints=500)

qemf_fw = Firework([qemft], name="QeMeanField")

wf = BgwWorkflow(qemf_fw, name="QE Bandstructure Task")

wf.to_file('gaas_bs_plot.yaml')

