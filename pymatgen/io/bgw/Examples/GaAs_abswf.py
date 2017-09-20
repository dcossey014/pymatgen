#!/usr/bin/env python

from fireworks import Firework, Workflow, LaunchPad

import os, copy, string, json, pprint
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.bgw.interfaces import QeMeanFieldTask, BgwFirework, BgwWorkflow
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
    commands.getoutput("/p/home/apps/ccm/opt/BerkeleyGW-1.1-beta2/Visual/gsphere.py gsphere.inp gsphere.out")
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

# defaults for this script
prev_qemf_dir=None

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

###
### below are model dependenant parameters that are required by BGW below
###
finder = SpacegroupAnalyzer(s)
s_prim = finder.get_primitive_standard_structure()
"""
Gives a structure with a primitive cell according to certain standards
the standards are defined in Setyawan, W., & Curtarolo, S. (2010).
High-throughput electronic band structure calculations:
Challenges and tools. Computational Materials Science,
49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010

Returns:
    The structure in a primitive standardized cell
"""
print s_prim
print s_prim.__repr__()

print "##########################################################"
print "### the following variables are set according to the model"
print "### and the parameters above."
print "###"
# get rid of spaces for pwscf and pw2bgw file prefixes
formula_prefix=s_prim.formula.replace(' ','')

# number_bands may later be attenuated to accomodate larger cuttoffs for more G vecs
# in sigma and epsilon
fftdum,number_bands=gsphere(screening_cutoff)

print "screening_cutoff =",screening_cutoff," number_bands =",number_bands
num_val_bands=int(num_val_electrons/2) # if this is not naturally an integer, then we should be doing spin polarized dft
print "num_val_electrons=",num_val_electrons
print "num_val_bands=",num_val_bands

krn_number_val_bands=num_val_bands-qp_band_index_min+1  
krn_number_cond_bands=qp_band_index_max - num_val_bands 

print "krn_number_val_bands =", krn_number_val_bands
print "krn_number_cond_bands =", krn_number_cond_bands

abs_number_val_bands_coarse = krn_number_val_bands
abs_number_cond_bands_coarse = krn_number_cond_bands

print "abs_number_val_bands_coarse =", abs_number_val_bands_coarse
print "abs_number_cond_bands_coarse =", abs_number_cond_bands_coarse

abs_number_val_bands_fine = abs_number_val_bands_coarse - n_val_bands_reduct
abs_number_cond_bands_fine = abs_number_cond_bands_coarse - n_cond_bands_reduct

print "abs_number_val_bands_fine =", abs_number_val_bands_fine
print "abs_number_cond_bands_fine =", abs_number_cond_bands_fine 

# choose number of bands to calculate based on above scheme with 
# a pad of two for degeneracy check
nbnd_wfn=number_bands+degeneracy_pad
nbnd_wfn_co=number_bands+degeneracy_pad 
nbnd_wfnq=num_val_bands+degeneracy_pad
nbnd_wfnq_fi=nbnd_wfnq
nbnd_wfn_fi=num_val_bands+abs_number_cond_bands_fine+degeneracy_pad

print "number of bands for WFN =", nbnd_wfn
print "number of bands for WFN_co =", nbnd_wfn_co
print "number of bands for WFNq =", nbnd_wfnq
print "number of bands for WFNq_fi =", nbnd_wfnq_fi
print "number of bands for WFN_fi =", nbnd_wfn_fi



###
### more workflow dependant setup that can be mostly hidden from this level
###

fftw_grid, dummy=gsphere(e_cut_mean_field)
print "fftw_grid from gsphere = ",fftw_grid
print ""

if cmplx_real == 'cmplx':
    cmplx_bool=True
else:
    cmplx_bool=False

kpoints = {}
kpoints['scf'] = kpts_co
kpoints['wfn'] = kpts_co
kpoints['wfn_co'] = kpts_co
kpoints['wfn_fi'] = kpts_fi
kpoints['wfnq'] = kpts_co
kpoints['wfnq_fi'] = kpts_fi


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
qemft_control['scf']=scf_control
for i in ['wfn', 'wfnq', 'wfn_co', 'wfn_fi', 'wfnq_fi']:
    qemft_control[i] = copy.deepcopy(wfn_control)

scf_system={'ecutwfc': e_cut_mean_field }

qemft_system = {}
for i in ['scf', 'wfn', 'wfnq', 'wfn_co', 'wfn_fi', 'wfnq_fi']:
    qemft_system[i] = copy.deepcopy(scf_system)

qemft_system['wfn']['nbnd'] = nbnd_wfn
qemft_system['wfnq']['nbnd'] = nbnd_wfnq
qemft_system['wfn_co']['nbnd'] = nbnd_wfn_co
qemft_system['wfn_fi']['nbnd'] = nbnd_wfn_fi
qemft_system['wfnq_fi']['nbnd'] = nbnd_wfnq_fi

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
qe_wfn_co_in_pw2bgw['vxc_diag_nmax']=number_bands
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmin']=0
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmax']=0

qemft = QeMeanFieldTask(structure=s_prim, kpoints=kpoints, pseudo_dir=pseudo_dir,
        mpi_cmd=mpi_cmd, pw_cmd=pw_cmd, pw2bgw_cmd=pw2bgw_cmd, qshift=qshift,
        fftw_grid=fftw_grid, qe_control=qemft_control, qe_system=qemft_system,
        qe_electrons=qemft_electrons, qe_pw2bgw=qemft_pw2bgw, mf_tasks=mean_field_tasks)

qemf_fw = Firework([qemft], name="QeMeanField")

eps_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
            kpoints=kpts_co, qshift=qshift, filename='epsilon.inp',
            qemf_dir=prev_qemf_dir)

eps_inp.epsilon_cutoff=screening_cutoff
eps_inp.number_bands = number_bands
namei="Epsilon Task"
eps_fw = BgwFirework(eps_inp, name=namei, complex=cmplx_bool, mpi_cmd=mpi_cmd)

#eps_fw.add_fw_to_launchpad()

sig_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
                kpoints=kpts_co, qshift=qshift, filename='sigma.inp',
                qemf_dir=prev_qemf_dir)

#sig_inp.screened_coulomb_cutoff = 10.0 # defaults to epsilon_cutoff
#sig_inp.bare_coulomb_cutoff = screening_cutoff # defaults to wfn cutoff
#sig_inp.band_index_min = 1

# must specify number_bands
sig_inp.number_bands=number_bands
sig_inp.screening_semiconductor = True
sig_inp.band_index_min=qp_band_index_min
sig_inp.band_index_max=qp_band_index_max

# ugh, sigh
eqp = '/app/ccm/opt/BerkeleyGW-1.2.0-hdf5/bin/eqp.py'
ppx = ' '.join([eqp, 'eqp1', './sigma_hp.log', './eqp_co.dat'])

sig_fw = BgwFirework(sig_inp, name="Sigma Task", ppx=ppx, complex=cmplx_bool, mpi_cmd=mpi_cmd)

bgw_wf=BgwWorkflow(qemf_fw, eps_fw, sig_fw, name="Sigma Converge")

krn_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
                qshift=qshift, filename='kernel.inp',
                qemf_dir=prev_qemf_dir)

# these are set above as part of "required" workflow parameters
krn_inp.screened_coulomb_cutoff = screening_cutoff
krn_inp.number_cond_bands = krn_number_cond_bands 
krn_inp.number_val_bands = krn_number_val_bands

# these are defaults for the type of run we are doing
krn_inp.use_symmetries_coarse_grid = True
krn_inp.screening_semiconductor = True

krn_fw = BgwFirework(krn_inp, name="Kernel Task", complex=cmplx_bool,  mpi_cmd=mpi_cmd)

abs_inp = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
                qshift=qshift, filename='absorption.inp',
                qemf_dir=prev_qemf_dir)

abs_inp.number_val_bands_fine = abs_number_val_bands_fine
abs_inp.number_val_bands_coarse = abs_number_val_bands_coarse
abs_inp.number_cond_bands_fine = abs_number_cond_bands_fine 
abs_inp.number_cond_bands_coarse = abs_number_cond_bands_coarse

# these are defaults for the type of run are doing
abs_inp.diagonalization = True
abs_inp.use_symmetries_coarse_grid = True
abs_inp.no_symmetries_fine_grid = True
abs_inp.no_symmetries_shifted_grid = True
abs_inp.screening_semiconductor = True
abs_inp.use_velocity = True
abs_inp.gaussian_broadening = True
abs_inp.energy_resolution = 0.15
abs_inp.eqp_co_corrections = True

abs_fw = BgwFirework(abs_inp, name="Absorption Task", complex=cmplx_bool,  mpi_cmd=mpi_cmd)

bgw_wf2 = BgwWorkflow(qemf_fw, eps_fw, sig_fw, krn_fw, abs_fw, 
                    name="QE/BGW Optical Task")

bgw_wf2.preserve_worker()
bgw_wf2.to_file('GaAs_abswf.yaml')
