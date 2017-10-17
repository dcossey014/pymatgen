#!/usr/bin/env python
from __future__ import division, unicode_literals

import unittest
import os

from pymatgen import Structure
from pymatgen.io.espresso.inputs import PWInput, PWInputError, PWOutput
from pymatgen.io.espresso.custodian_jobs import WritePwscfInputTask, SimplePWTask
from fireworks import Firework, Workflow, LaunchPad

test_dir="/home/gkedz/python/MaterialsProject/pymatgen/test_files"

s=Structure.from_file('./mp-149_Si.cif')
my_control={"calculation": "scf", "pseudo_dir": "/home/gkedz/BGW/WorkFlows/PWSCF"}
my_pseudo={"Si": "Si.UPF"}
my_system={"ecutwfc": 25}
#my_kpoints_grid=(5,5,5)

with open('kpoints','r') as kin:
    my_kpoints=kin.readlines()
kpoints_mode=my_kpoints.pop(0).split()[-1]

pw=PWInput(s,control=my_control,pseudo=my_pseudo,system=my_system,
                kpoints_grid=my_kpoints,kpoints_mode=kpoints_mode)


pw.write_file("si_pw.in")

