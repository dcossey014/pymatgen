#!/usr/bin/env python

from fireworks import Firework, Workflow, LaunchPad
from pymatgen.io.espresso.inputs import PwInput

import os, copy, string, json, pprint
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.bgw.interfaces import BgwFirework, BgwWorkflow
from pymatgen.io.espresso.interfaces import QeMeanFieldTask
from pymatgen.io.bgw.inputs import BgwInput
import numpy as np
import commands
import re

# define gsphere here temporarily, move to appropriate place later, 
# but I think ecutwfn in pw.x and cuttoffs and bands in epsilon and sigma
# need to be consistent.  Also what about the rho cuttoff?
def gsphere(e_cut):
    #print "in gsphere e_cut=",e_cut
    gsin=open("gsphere.inp",'w')
    a2b=1.88973
    cell_vecs=np.asarray(s_prim.lattice.matrix,dtype=float)
    cell_vecs_bohr=cell_vecs*a2b
    for vec in cell_vecs_bohr:
        gsin.write( " %14.8f %14.8f %14.8f\n" % (vec[0], vec[1], vec[2]) )
    gsin.write("%5.1f\n" % e_cut)
    gsin.write("0 0 0\n")
    gsin.write("false\n")
    gsin.write("0 0 0\n")
    gsin.close()
    commands.getoutput("cat gsphere.inp")
    #commands.getoutput("/p/home/apps/ccm/opt/BerkeleyGW-1.1-beta2/Visual/gsphere.py gsphere.inp gsphere.out")
    commands.getoutput("/Users/pettt/tmp/Visual/gsphere.py gsphere.inp gsphere.out")
    gsout=open("gsphere.out",'r')
    line=gsout.readline()
    while line:
        line=gsout.readline()
        if "grid =" in line:
            line=gsout.readline()
            line=gsout.readline()
            #print line
            match=re.search(r"\(\s+(\d+)\s+(\d*)\s+(\d*)\s+\) -- factors\s+(.+)",line)
            if match:
                fft_x=int(match.group(1))
                fft_y=int(match.group(2))
                fft_z=int(match.group(3))
                factstr=match.group(4)
                #print "fft_x=",fft_x,"fft_y=",fft_y,"fft_z=",fft_z
                #print "factstr=",factstr
                factors=factstr.split(',')
                #print "factors[0]=",factors[0]
            break
    
    while line:
        line=gsout.readline()
        #print line
        if "ng =" in line:
            match=re.search(r"ng =\s+(\d+)",line)
            if match:
                number_bands=int(match.group(1))
                #print "number_bands=",number_bands
            break
    
    grid_accuracy=int(factors[0])  
    fftw_grid = [fft_x*grid_accuracy,fft_y*grid_accuracy,fft_z*grid_accuracy]
    #print "fftw_grid=",fftw_grid
    gsout.close()
    return fftw_grid, number_bands

###
### USER PARAMETERS. These are the most important parameters to examine carefully 
### for this particular model and workflow.
###

s = Structure.from_file('GaAs_sg216_exp.cif')

# This is adjusted according to the number of cpus used
pw_cmd = '/app/espresso/platforms/espresso-6.1.0/bin/pw.x -npool 4 -ndiag 36'
pw2bgw_cmd = '/app/espresso/platforms/espresso-6.1.0/bin/pw2bgw.x -npool 4 -ndiag 36'
#prev_qemf_dir='/workspace/gkedz/GaAs/old_blk_20170818/launcher_2017-08-07-18-31-20-782413/ESPRESSO'
mpi_cmd = 'mpiexec_mpt -n 72'
pseudo_dir = '/app/ccm/PPs/ONCVPSP/QE'
kpts_co = [6,6,6]
kpts_fi= [8,8,8]
e_cut_mean_field=60.0 
screening_cutoff=10.0 # screening energy cutoff 
cmplx_real='cmplx'

# this is the quasi-particle band range for the self-energy and dispersion matrices
qp_band_index_min=6   # we are not interested in inner d states for now
qp_band_index_max=17  # this should be consistent with nbands in wfnq and wfnq_fi

# these are the reduction terms for the fine grid in absorption
n_val_bands_reduct=2  # do not know guidance here yet, just guess
n_cond_bands_reduct=2 # do not know guidance here yet, just guess

# these will be set automatically or by default laters

# this helps get desired number of bands in WFN* after degeneracy check
degeneracy_pad=2
num_val_electrons=18 # get this from pseudo potential files
qshift=[0.0, 0.0, 0.001]

print "##########################################################"
print "### User parameters for GW BSE calculation"
print "###"
print ""
print "kpts_co =", kpts_co
print "kpts_fi =", kpts_fi
print "e_cut_mean_field =", e_cut_mean_field
print "screening_cutoff =",screening_cutoff
print "state representation field type cmplx_real =",cmplx_real
print "qp_band_index_min =",qp_band_index_min
print "qp_band_index_max =",qp_band_index_max
print "n_val_bands_reduct =",n_val_bands_reduct
print "n_cond_bands_reduct =",n_cond_bands_reduct
print ""

# Needed only for gsphere
finder = SpacegroupAnalyzer(s)
s_prim = finder.get_primitive_standard_structure()

# number_bands may later be attenuated to accomodate larger cuttoffs for more G vecs
# in sigma and epsilon
fftdum,number_bands=gsphere(screening_cutoff)

# Three Different Cases for config_file: Absolute Path, 
# Relative Path, and Automatic Search in HOME directory
config_file = '/Users/pettt/tmp/pymatgen/pymatgen/io/espresso/espresso_interface_defaults.yaml'
#config_file = 'espresso_interface_defaults2.yaml'
#config_file = None


# 2 Cases of Kpoints: Automatic and Crystal
# Automatic Case
kps_mode = 'automatic'
kps = [8,8,8]
kps_shift = [0.5, 0.5, 0.5]
qshift = [0,0,0]

# Crystal Case
"""
kps_mode = 'crystal'
kps = ['  120',
        '  0.062500000  0.062500000  0.062500000   1.0\n',
        '  0.062500000  0.062500000  0.187500000   3.0\n',
        '  0.062500000  0.062500000  0.312500000   3.0\n',
        '  0.062500000  0.062500000  0.437500000   3.0\n'
        ]
kps_shift = [0,0,0]
"""

qe_scf = PwInput(s, kpoints_mode=kps_mode, kpoints_grid=kps, 
                kpoints_shift=kps_shift, config_file=config_file)
qe_scf.system['ecutwfc'] = 30

qe_scf.write_file('scf.in')

