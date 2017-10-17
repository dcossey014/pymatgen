#!/usr/bin/env python

from fireworks import Firework, Workflow, LaunchPad

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


### begin user input such as it is.  It should be one command
s = Structure.from_file('GaAs_sg216_exp.cif')
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

    

# get rid of spaces for pwscf and pw2bgw file prefixes
formula_prefix=s_prim.formula.replace(' ','')

### user parameters
# This is adjusted according to the number of cpus used
pw_cmd = '/app/espresso/platforms/espresso-6.1.0/bin/pw.x -npool 4 -ndiag 36'
pw2bgw_cmd = '/app/espresso/platforms/espresso-6.1.0/bin/pw2bgw.x -npool 4 -ndiag 36'
#prev_qemf_dir='/workspace/gkedz/GaAs/block_2017-07-29-12-41-31-296723/launcher_2017-07-29-13-18-14-601144/ESPRESSO'
mpi_cmd = 'mpiexec_mpt -n 144'
pseudo_dir = '/app/ccm/PPs/ONCVPSP/QE'
e_cut_mean_field=200.0
kpts_co = [6,6,6]
kpts_fi= [8,8,8]
qshift=[0.0, 0.0, 0.001]
cmplx_real='cmplx'
if cmplx_real == 'cmplx':
    cmplx_bool=True
else:
    cmplx_bool=False

### workflow dependant setup
fftw_grid, dummy=gsphere(e_cut_mean_field)
print "fftw_grid from gsphere = ",fftw_grid
eps_inp=[None]*6
eps_fw=[None]*6
nbnds=[]
i=0
cut_min=5
cut_max_plus=35
cut_incr=5
for icut in range(cut_min,cut_max_plus,cut_incr):
    cut=float(icut)
    fftdum,number_bands=gsphere(cut)
    print "cut =",cut," number_bands =",number_bands

    eps_inp[i] = BgwInput(s_prim, pseudo_dir=pseudo_dir, cmplx_real=cmplx_real,
                kpoints=kpts_co, qshift=qshift, filename='epsilon.inp')

    eps_inp[i].epsilon_cutoff=cut
    eps_inp[i].number_bands = number_bands
    nbnds.append(number_bands)
    namei="Epsilon Task {}".format(i)
    eps_fw[i] = BgwFirework(eps_inp[i], name=namei, complex=cmplx_bool)
    i=i+1

num_eps_fw=i
print "num_eps_fw=",num_eps_fw

nbnd_wfn=max(nbnds)
nbnd_wfnq=16
nbnd_wfn_co=45
nbnd_wfn_fi=16
nbnd_wfnq_fi=16

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
qe_wfn_co_in_pw2bgw['vxc_diag_nmax']=45
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmin']=0
qe_wfn_co_in_pw2bgw['vxc_offdiag_nmax']=0

qemft = QeMeanFieldTask(structure=s_prim, kpoints=kpoints, pseudo_dir=pseudo_dir,
        mpi_cmd=mpi_cmd, pw_cmd=pw_cmd, pw2bgw_cmd=pw2bgw_cmd, qshift=qshift,
        fftw_grid=fftw_grid, qe_control=qemft_control, qe_system=qemft_system,
        qe_electrons=qemft_electrons, qe_pw2bgw=qemft_pw2bgw, mf_tasks=mean_field_tasks)

qemf_fw = Firework([qemft], name="QeMeanField")

bgw_wf = BgwWorkflow(qemf_fw, deps_dict={qemf_fw: []}, name="QE/Converge Epsilon")

for i in range(0,num_eps_fw):
    namei="Epsilon Task {}".format(i)
    eps_cnv_deps={ qemf_fw: [eps_fw[i]] }
    bgw_wf.add_fw(eps_fw[i], deps=eps_cnv_deps)
    i=i+1

bgw_wf.preserve_worker()
bgw_wf.to_file('gaas_cnv_loop3.yaml')
