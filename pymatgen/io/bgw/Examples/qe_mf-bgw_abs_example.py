#!/usr/bin/env python

from fireworks import Firework, Workflow

import os
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.bgw.interfaces import QeMeanFieldTask, BgwAbsTask

s = Structure.from_file('/app/ccm/Si_symm.cif')
finder = SpacegroupAnalyzer(s)

s_prim = finder.get_primitive_standard_structure()
s_prim.translate_sites([0,1], [-0.125,-0.125,-0.125], to_unit_cell=False)
print s_prim
print s_prim.__repr__()


kpoints = [5,5,5]
qshift=[0.0, 0.0, 0.001]
fftw_grid = [24,24,24]
os.environ['KGRID_EXEC'] = '/apps/ccm/opt/BerkeleyGW-1.1-beta2/bin/kgrid.x'

pseudo_dir = '/app/ccm/PPs/ONCVPSP/QE'
mpi_cmd = 'mpiexec_mpt -n 36'

pw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw.x'
pw2bgw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw2bgw.x'

qemft_control = {'scf': {"calculation" : "scf", "pseudo_dir": pseudo_dir, 
                    'prefix': s_prim.formula.replace(' ', ''),
                    'tstress': True, 'tprnfor': True},
                'wfn': {"calculation" : "bands", "pseudo_dir": pseudo_dir,
                    'prefix': s_prim.formula.replace(' ', '')},
                'wfnq': {"calculation" : "bands", "pseudo_dir": pseudo_dir,
                    'prefix': s_prim.formula.replace(' ', '')},
                'wfn_co': {"calculation" : "bands", "pseudo_dir": pseudo_dir,
                    'prefix': s_prim.formula.replace(' ', '')},
                'wfn_fi': {"calculation" : "bands", "pseudo_dir": pseudo_dir,
                    'prefix': s_prim.formula.replace(' ', '')},
                'wfnq_fi': {"calculation" : "bands", "pseudo_dir": pseudo_dir,
                    'prefix': s_prim.formula.replace(' ', '')}
                }

qemft_system = {'scf': {'ecutwfc': 25.0 },
                'wfn': {'ecutwfc': 25.0 , 'nbnd': 30},
                'wfnq': {'ecutwfc': 25.0    , 'nbnd': 7},
                'wfn_co': {'ecutwfc': 25.0  , 'nbnd': 33},
                'wfn_fi': {'ecutwfc': 25.0  , 'nbnd': 11},
                'wfnq_fi': {'ecutwfc': 25.0 , 'nbnd': 4}
                }

qemft_electrons = {'scf': {'conv_thr': 1.0e-10, 'diago_full_acc': True},
                'wfn': {'conv_thr': 1.0e-10, 'diago_full_acc': True, 'startingwfc': 'random'},
                'wfnq': {'conv_thr': 1.0e-10, 'diago_full_acc': True, 'startingwfc': 'random'},
                'wfn_co': {'conv_thr': 1.0e-10, 'diago_full_acc': True, 'startingwfc': 'random'},
                'wfn_fi': {'conv_thr': 1.0e-10, 'diago_full_acc': True, 'startingwfc': 'random'},
                'wfnq_fi': {'conv_thr': 1.0e-10, 'diago_full_acc': True, 'startingwfc': 'random'}
                }

qemft = QeMeanFieldTask(structure=s_prim, kpoints=kpoints, pseudo_dir=pseudo_dir,
        mpi_cmd=mpi_cmd, pw_cmd=pw_cmd, pw2bgw_cmd=pw2bgw_cmd, qshift=qshift,
        fftw_grid=fftw_grid, qe_control=qemft_control, qe_system=qemft_system,
        qe_electrons=qemft_electrons)

qemft.run_task(fw_spec={})

fw1 = Firework([qemft], name="QeMeanFieldExample")
#fw1.to_file("qemft_fw_example.yaml")

#Begin BGW Calculation
bgw_dir = "/apps/ccm/opt/BerkeleyGW-1.1-beta2/bin/"
eps = bgw_dir+'epsilon.real.x'
sig = bgw_dir+'sigma.real.x'
eqp = bgw_dir+'eqp.py'
krn = bgw_dir+'kernel.real.x'
absorp = bgw_dir+'absorption.real.x'

bgw_params = {'epsilon': {'epsilon_cutoff': 10.0, 'number_bands': 29},
            'sigma': {'screened_coulomb_cutoff': 10.0, 'bare_coulomb_cutoff': 25.0,
                    'number_bands': 29, 'band_index_min': 1, 'band_index_max':14,
                    'screening_semiconductor': ''},
            'kernel': {'number_cond_bands': 10, 'screened_coulomb_cutoff': 10.0,
                    'bare_coulomb_cutoff': 25.0, 'use_symmetries_coarse_grid': '',
                    'screening_semiconductor': ''},
            'absorption':{'diagonalization': '', 'number_cond_bands_coarse': 10,
                    'number_cond_bands_fine': 6, 'use_symmetries_coarse_grid': '',
                    'no_symmetries_fine_grid': '', 'no_symmetries_shifted_grid': '',
                    'screening_semiconductor': '', 'use_velocity': '',
                    'gaussian_broadening': '', 'energy_resolution': 0.15,
                    'eqp_co_corrections': ''}}

bgwabs = BgwAbsTask(structure=s_prim, mpi_cmd=mpi_cmd, eps_exe=eps, sig_exe=sig,
                eqp_exe=eqp, krn_exe=krn, abs_exe=absorp, bgw_isp=bgw_params,
                kpoints=kpoints, pseudo_dir=pseudo_dir, qshift=qshift)

bgwabs.run_task(fw_spec={})

fw2 = Firework([bgwabs], name='BGW Absorption')
fw2.to_file('bgw_fw_example.yaml')

