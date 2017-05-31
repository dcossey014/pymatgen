#!/usr/bin/env python

from fireworks import Firework

import os
from pymatgen import Structure
from pymatgen.io.bgw.interfaces import QeMeanFieldTask

s = Structure.from_file(os.environ['HOME']+'/mp-149_Si.cif')
s_prim = s.get_primitive_structure()

kpoints = [5,5,5]
fftw_grid = [24,24,24]
os.environ['KGRID_EXEC'] = '/apps/ccm/opt/BerkeleyGW-1.1-beta2/bin/kgrid.x'

pseudo_dir = '/app/ccm/PPs/ONCVPSP/QE'
mpi_cmd = 'mpiexec_mpt -n 36'

pw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw.x'
pw2bgw_cmd = '/app/espresso/platforms/espresso-5.3.0/bin/pw2bgw.x'

qemft = QeMeanFieldTask(structure=s_prim, kpoints=kpoints, pseudo_dir=pseudo_dir,
        mpi_cmd=mpi_cmd, pw_cmd=pw_cmd, pw2bgw_cmd=pw2bgw_cmd, fftw_grid=fftw_grid)

qemft.run_task(fw_spec={})

fw1 = Firework([qemft], name="QeMeanFieldExample")
fw1.to_file("qemft_fw_example.yaml")
