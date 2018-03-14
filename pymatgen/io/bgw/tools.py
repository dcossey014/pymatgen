__author__ = "David Cossey"
__copyright__ = "n/a"
__version__ = "0.0.1a"
__maintainer__ = "David Cossey"
__email__ = "dcossey014@gmail.com"
__date__ = "2018-01-31"

import six, glob, sys, errno
import os, abc, commands
import mmap, fnmatch, re, subprocess
import copy as cp
from bisect import bisect_left
import numpy as np

from pymatgen.io.bgw.kgrid import Kgrid, generate_kpath

from monty.io import zopen
from monty.dev import deprecated
from monty.json import MSONable
from monty.serialization import loadfn
from collections import defaultdict, OrderedDict
from fireworks import FireTaskBase, explicit_serialize, FWAction
from pymatgen import Structure


def gsphere(structure, e_cut):
    s = structure if isinstance(structure, Structure) else Structure.from_dict(
                structure)

    #print "in gsphere e_cut=",e_cut
    gsin=open("gsphere.inp",'w')
    a2b=1.88973
    cell_vecs=np.asarray(s.lattice.matrix,dtype=float)
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
    #commands.getoutput("/Users/pettt/tmp/Visual/gsphere.py gsphere.inp gsphere.out")
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

def get_num_val_electrons(struct,pseudo_dir=None,config_file=None):

    def either(c):
        return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c

    if not pseudo_dir:
        if config_file:
            config_dict=loadfn(config_file)
        elif os.path.exists(os.path.join(os.environ['HOME'],'bgw_interface_defaults.yaml')):
            config_dict = loadfn(os.path.join(os.environ['HOME'],
            'bgw_interface_defaults.yaml'))
            pseudo_dir=config_dict['espresso']['control']['pseudo_dir']
        else:
            print "Don't have a pseudpotential directory or config file"
            exit()

    if not pseudo_dir:
        print "Can not find pseudo potential files\n"

    comp=struct.composition
    elem_amt_dict=comp.get_el_amt_dict()

    num_val_e=0
    for i in struct.symbol_set:
        amt=elem_amt_dict.get(str(i))
        pattern="{}/{}_*.UPF".format(pseudo_dir,i)
        new_pattern=''.join(either(char) for char in pattern)
        pps=glob.glob(new_pattern)
        ppfile=pps[-1]
        ppfh=open(ppfile)
        for line in ppfh:
            if "z_valence" in line:
                l=line.strip().split()
                z_val=float(l[-1][:-1])
                break
        num_val_e=num_val_e+z_val*amt

    return int(num_val_e)
