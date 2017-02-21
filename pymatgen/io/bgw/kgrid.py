import os
import sys
import string
import subprocess
import seekpath

import numpy as np
from enum import Enum
from pymatgen import Structure
from monty.serialization import loadfn

def generate_kpath(s, npts):
    pos = s.frac_coords
    nums = s.atomic_numbers
    cell = s.lattice.matrix
    struc = (cell, pos, nums)
    kpath = seekpath.get_path(struc)
    pcoords = kpath['point_coords']
    path = list(kpath['path'])
    rec_latt = kpath['reciprocal_primitive_lattice']

    #Generate a continuous path
    ep=0
    for i,p in enumerate(path):
        if ep != p[0] and i != 0:
            path.insert(i, (ep, p[0]))
        sp, ep = path[i][0], path[i][1]

    # Gather total distance of path to subdivide
    magnitude = []
    for i in path:
        vec = np.array(pcoords[i[1]]) - np.array(pcoords[i[0]])
        abc_vec = vec.dot(rec_latt)
        magnitude.append(np.linalg.norm(abc_vec))
    total_length = sum(magnitude)

    # Number of points per segment is divided evenly along path
    # Path is stored in gen_kpath and returned
    gen_kpath = [pcoords[path[0][0]]]
    for i,p in enumerate(path):
        length = np.array(pcoords[p[1]]) - np.array(pcoords[p[0]])
        num_divs = int(round(magnitude[i] / total_length * npts, 0))
        sub_div = length / num_divs
        for j in range(num_divs):
            gen_kpath.append(gen_kpath[-1] + sub_div)
    return gen_kpath


class QeMeanFieldGrids(object):
    '''
    Docstring
    '''
    def __init__(self, structure, kpoints=[5,5,5],
            offset_type="Monkhorst-Pack",
            qshift=[0.000, 0.000, 0.001],
            fftw_grid=[0,0,0], bgw_rev_off=False,
            log_cart_kpts=False):

        self.structure = structure
        self.kpoints = ( kpoints if isinstance(kpoints, dict) \
                            else {'scf': kpoints, 'wfn': kpoints, 
                                'wfn_co': kpoints, 'wfnq': kpoints})
        self.offset_type = ( offset_type if isinstance(offset_type, dict) \
                            else {'scf': "Monkhorst-Pack", 
                                'wfn': "Monkhorst-Pack", 'wfn_co': 'Gamma', 
                                'wfnq': "Monkhorst-Pack", 
                                'wfn_fi': [0.47, 0.37, 0.31],
                                'wfnq_fi': [0.47, 0.37, 0.32]} )
        self.qshift = qshift
        self.fftw_grid = fftw_grid
        self.bgw_rev_off = bgw_rev_off
        self.log_cart_kpts = log_cart_kpts
        self.grids={}
        
        if 'wfn_co' not in self.kpoints.keys():
            self.kpoints['wfn_co'] = ( [ i for i in \
                        self.kpoints['wfn'] ] )
        elif self.kpoints['wfn_co'] != self.kpoints['wfn']:
            print("WFN_co and WFN kpoints must be equal")
            print("Resetting WFN_co kpoints to WFN kpoints")
            self.kpoints['wfn_co'] = [ i for i in self.kpoints['wfn'] ]

        if 'wfn_fi' not in self.kpoints.keys(): 
            self.kpoints['wfn_fi'] = ( [ int(i*2) for i in \
                        self.kpoints['wfn'] ] )

        if 'wfnq' not in self.kpoints.keys():
            self.kpoints['wfnq'] = self.kpoints['wfn']
        elif self.kpoints['wfn'] != self.kpoints['wfnq']:
            print("WFNq and WFN kpoints must be equal")
            print("Resetting WFNq kpoints to WFN kpoints")
            self.kpoints['wfnq'] = [ i for i in self.kpoints['wfn'] ]

        if 'wfnq_fi' not in self.kpoints.keys():
            self.kpoints['wfnq_fi'] = self.kpoints['wfn_fi']
        elif self.kpoints['wfnq_fi'] != self.kpoints['wfn_fi']:
            print("WFNq_fi and WFN_fi kpoints must be equal")
            print("Resetting WFNq_fi kpoints to WFN_fi kpoints")
            self.kpoints['wfnq_fi'] = [ i for i in self.kpoints['wfn_fi'] ]


    def generate_kgrid(self,qe_task='scf'):
        #print("running task: {}".format(qe_task))

        #case scf
        if qe_task == 'scf':
            kpoints=self.kpoints['scf']
            offset_type=self.offset_type['scf']
            qshift=[0.0,0.0,0.0]

        #case wfn
        elif qe_task == 'wfn':
            kpoints=self.kpoints['wfn']
            offset_type=self.offset_type['wfn']
            qshift=[0.0,0.0,0.0]

        #case wfn_co
        elif qe_task == 'wfn_co':
            kpoints=self.kpoints['wfn_co']
            offset_type=self.offset_type['wfn_co']
            qshift=[0.0,0.0,0.0]           
        
        #case wfn_fi
        elif qe_task == 'wfn_fi':
            kpoints=self.kpoints['wfn_fi']
            offset_type=self.offset_type['wfn_fi']
            qshift=[0.0,0.0,0.0]

        #case wfnq
        elif qe_task == 'wfnq':
            kpoints=self.kpoints['wfnq']
            offset_type=self.offset_type['wfnq']
            qshift=self.qshift

        #case wfnq_fi
        elif qe_task == 'wfnq_fi':
            kpoints=self.kpoints['wfnq_fi']
            offset_type=self.offset_type['wfnq_fi']
            qshift=[0.0,0.0,0.0]
        else:
            print "Unknown QE task in QeMeanFieldGrids"
            exit()

        self.grids[qe_task] = Kgrid(self.structure, kpoints=kpoints, offset_type=offset_type, 
                qshift=qshift, fftw_grid=self.fftw_grid, bgw_rev_off=self.bgw_rev_off, 
                log_cart_kpts=self.log_cart_kpts)

        return self.grids[qe_task].generate_kpoints(qe_task)


class Kgrid(object):
    '''
    Docstring
    '''
    def __init__(self, structure, kpoints=[1,1,1],
            offset_type="Monkhorst-Pack", 
            qshift=[0.000, 0.000, 0.000], 
            fftw_grid=[0, 0, 0], bgw_rev_off=False,
            log_cart_kpts=False):

        self.structure = structure
        self.kpoints = kpoints
        self.offset_type = offset_type
        self.qshift = qshift
        self.fftw_grid = fftw_grid
        self.bgw_rev_off = bgw_rev_off
        self.log_cart_kpts = log_cart_kpts

    def write_input(self, filename):
        prim_struct = self.structure.get_primitive_structure()
        prim_lvs = prim_struct.lattice.matrix
        prim_coords = prim_struct.cart_coords
        if isinstance(self.offset_type, list):
            kgrid_offset = "{0:< 5.3f} {1:< 5.3f} {2:< 5.3f}\n".format(
                    self.offset_type[0], self.offset_type[1], 
                    self.offset_type[2])
        elif self.offset_type.lower() == 'gamma':
            kgrid_offset = "{0:< 5.3f} {0:< 5.3f} {0:< 5.3f}\n".format(0.0)
        elif self.offset_type.lower() == "random":
            kgrid_offset = "{0:< 5.3f} {1:< 5.3f} {2:< 5.3f}\n".format(0.47, 
                                                            0.37, 0.32)
        else:
            kgrid_offset = "{0:< 5.3f} {0:< 5.3f} {0:< 5.3f}\n".format(0.5)

        element = Enum('Elem', ",".join(prim_struct.symbol_set).encode('ascii', 'ignore'))

        with open(filename, 'w') as fout:
            fout.write("{:< 5} {:< 5} {:< 5}\n".format(self.kpoints[0], 
                                    self.kpoints[1], self.kpoints[2]))
            fout.write("{}".format(kgrid_offset))
            fout.write("{:< 5} {:< 5} {:< 5}\n\n".format(self.qshift[0], 
                                    self.qshift[1], self.qshift[2]))

            for lv in prim_lvs:
                fout.write("{:< 12.8f} {:< 12.8f} {:< 12.8f}\n".format(lv[0], 
                                                            lv[1], lv[2]))
            fout.write("{}\n".format(prim_struct.num_sites))

            for i,el in enumerate(prim_struct.species):
                fout.write("{} ".format(element[el.__str__()].value))
                crds = prim_struct.cart_coords[i]
                fout.write("{:< 12.8f} {:< 12.8f} {:< 12.8f}\n".format(crds[0], 
                                                        crds[1], crds[2]))

            fout.write("{:< 5} {:< 5} {:< 5}\n".format(self.fftw_grid[0], 
                                    self.fftw_grid[1], self.fftw_grid[2]))

            fout.write("{}\n".format(".true." if self.bgw_rev_off else ".false."))
            fout.write("{}\n".format(".true." if self.log_cart_kpts else ".false."))

    def generate_kpoints(self, basename, config_file=None, kgridx=None):
        if config_file:
            config_dict = loadfn(config_file)
            kgrid_exec = os.path.join(config_dict['BGW_DIR'], 'kgrid.x')
        elif os.path.exists(os.path.join(os.environ['HOME'],
                                        'bgw_interface_defaults.yaml')):
            config_dict = loadfn(os.path.join(os.environ['HOME'],
                                        'bgw_interface_defaults.yaml'))
            kgrid_exec = os.path.join(config_dict['BGW_DIR'], 'kgrid.x')
        else:
            config_dict = {}


        if kgridx:
            kgrid_exec = kgridx
        elif os.environ.get('KGRID_EXEC', None):
            kgrid_exec = os.environ['KGRID_EXEC']

        try:
            if not kgridx and not config_file and not os.environ['KGRID_EXEC']:
                return "No location for kgrid.x was given.  Specify location in config_file ",\
                        "as kgridx or in method call with kgridx=</path/to/kgrid.x>"
            elif config_file and kgridx:
                if config_dict['kgridx']:
                    return "Please specify kgridx in the configuration file or in the method call,"\
                            "not both.  Found multiple locations for kgrid.x"
            else:
                print("Continuing with kgrid.x from {}\n".format(kgrid_exec))
        except:
            print("Continuing with kgrid.x from {}\n".format(kgrid_exec))

        self.write_input(basename+".in")
        p = subprocess.Popen([kgrid_exec, basename+'.in', basename+'.out', basename+'.log'])
        p.wait()

        with open(basename+'.out', 'r') as fin:
            kpoint_grid = fin.readlines()
        rm_line = kpoint_grid.pop(0)
        kpoint_grid[-1] = kpoint_grid[-1].rstrip()

        return kpoint_grid
