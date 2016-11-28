#!/usr/bin/env python
from __future__ import division, unicode_literals

import unittest
import os

from pymatgen import Structure
from pymatgen.io.pwscf import PWInput, PWInputError, PWOutput
from pymatgen.io.bgw.pwscf_tasks import WritePwscfInputTask, SimplePWTask
from fireworks import Firework, Workflow, LaunchPad

test_dir="/home/gkedz/python/MaterialsProject/pymatgen/test_files"

s=Structure.from_file('./Si_sg227_exp.cif')
my_control={"calculation": "scf", "pseudo_dir": "/home/gkedz/BGW/WorkFlows/PWSCF"}
my_pseudo={"Si": "Si.UPF"}
my_system={"ecutwfc": 25}
my_kpoints_grid=(5,5,5)

wpwit=WritePwscfInputTask(structure=s, pseudo=my_pseudo, control=my_control,
    system=my_system,kpoints_grid=my_kpoints_grid)
runpw=SimplePWTask()

fw=Firework([wpwit,runpw])

fw.to_file("test_pw.yaml")
