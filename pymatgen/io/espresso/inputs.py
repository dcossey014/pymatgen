# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module implements input and output processing from PWSCF.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "3/27/15"

import pprint
import six, math, abc, glob, sys, errno
import os, fnmatch, re, subprocess
from bisect import bisect_left

from pymatgen.io.bgw.kgrid import Kgrid, generate_kpath

from monty.io import zopen
from monty.dev import deprecated
from monty.json import MSONable
from monty.serialization import loadfn

from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from copy import deepcopy as dcopy
from monty.re import regrep
from collections import defaultdict, OrderedDict


class PwInput(MSONable):
    '''
    Base input file class for Quantum Espresso pw.x
    '''

    #TODO: Write __setattr__ function to make things nice in a dictionary
    #TODO: Make all __dict__ calls attributes.  No reason for __dict__ since
    #       this is not a FireTaskBase Class.
    def __init__(self, structure, pseudo_dir=None, control={'calculation': 'scf'}, 
                system={}, electrons={}, ions={}, cell={}, kpoints_mode='automatic', 
                kpoints_grid=[1,1,1], kpoints_shift=[1,1,1], 
                reduce_structure=False, config_file=None):
        """
        Initializes a PWSCF input file for using pw.x in QuantumEspresso.

        Args:
            structure (Structure):  Input Structure
            pseudo_dir (str):  String of the directory containing pseudopotential 
                    files.
            control (dict):  Control parameters. Refer to official PWSCF doc
                    on supported parameters. Default to {"calculation": "scf"}
            system (dict):  System parameters. Refer to official PWSCF doc
                    on supported parameters. Default to None, ( {} ).
            electrons (dict):  Electron parameters. Refer to official PWSCF doc
                    on supported parameters. Default to None, ( {} ).
            ions (dict):  Ions parameters. Refer to official PWSCF doc
                    on supported parameters. Default to None, ( {} ).
            cell (dict):  Cell parameters. Refer to official PWSCF doc
                    on supported parameters. Default to None, ( {} ).
            kpoints_mode (str):  Kpoints generation mode. Default to 'automatic'.
            kpoints_grid (list):  The kpoint grid. Default to [1, 1, 1].
            kpoints_shift (list):  The shift for the kpoints. Defaults to
                    [1,1,1].
            reduce_structure (bool):  If True, structure will be converted to 
                    its reduced primitive structure.
            config_file (str):  Configuration file from which to pull user default
                    values for control, system, electrons, and ions.  Default is 
                    to look for specified file in User's HOME directory but will 
                    also accept absolute paths for alternative locations.
        """

        # Check for Configuration Defaults File and apply those changes.
        #print("In PwInput: config_file: {}".format(config_file))
        self.__dict__['config_file'] = config_file
        if config_file and os.path.isabs(config_file):
            config_dict = loadfn(config_file)
        elif config_file:
            fpath = os.path.join(os.environ['HOME'], config_file)
            config_dict = loadfn(fpath)
        elif os.path.exists(os.path.join(os.environ['HOME'],
                    'espresso_interface_defaults.yaml')):
            fpath = os.path.join(os.environ['HOME'], 'espresso_interface_defaults.yaml')
            config_dict = loadfn(fpath)
        else:
            config_dict = {}

        (self.__dict__['control'],
            self.__dict__['system'],
            self.__dict__['electrons'],
            self.__dict__['ions'],
            self.__dict__['cell']) = {'calculation': 'scf'}, {}, {}, {}, {}

        # Set Defaults from Configuration File if one is found
        if config_dict:
            self.__dict__['control'].update(
                    config_dict.get('control', {}) )
            self.__dict__['system'].update(
                    config_dict.get('system', {}) )
            self.__dict__['electrons'].update(
                    config_dict.get('electrons', {}) )
            self.__dict__['ions'].update(
                    config_dict.get('ions', {}) )
            self.__dict__['cell'].update(
                    config_dict.get('cell', {}) )

        #print("values after config_file: {}".format(self.control))
        
        # Override Defaults with User Supplied Values
        self.__dict__['structure'] = structure
        self.__dict__['pseudo_dir'] = pseudo_dir
        self.__dict__['control'].update(control if control else {})
        self.__dict__['system'].update(system if system else {})
        self.__dict__['electrons'].update(electrons if electrons else {})
        self.__dict__['ions'].update(ions if ions else {})
        self.__dict__['cell'].update(cell if cell else {})
        self.__dict__['kpoints_mode'] = kpoints_mode.lower() 
        self.__dict__['kpoints_grid'] = kpoints_grid
        self.__dict__['kpoints_shift'] = kpoints_shift
        self.__dict__['reduce_structure'] = reduce_structure


        # If using reduced primite Structure, set that before moving forward
        if self.reduce_structure:
            self.__dict__['structure'] = self.convert2primitive(self.structure)

        # Set prefix name automatically if not already set in control dictionary
        formula_prefix = self.structure.composition.reduced_formula
        self.control['prefix'] = ( formula_prefix 
                            if not 'prefix' in self.control.keys() 
                            else self.control['prefix'] )
        
        # Overwrite PseudoPotential Directory from Deafults File if User supplied
        # one during the initialization of the class.
        if pseudo_dir:
            self.control['pseudo_dir'] = pseudo_dir

        # Set dictionary for Atom types to use particular pseudopotential files.
        # might want to have some user override on this maybe?
        #print("pseudo_dir: {}".format(self.control['pseudo_dir']))
        def either(c):
            return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
        self.__dict__['pseudo'] = {}
        for i in self.structure.symbol_set:
            pattern="{}/{}_*.UPF".format(self.control['pseudo_dir'],i)
            new_pattern=''.join(either(char) for char in pattern)
            pps=glob.glob(new_pattern)
            self.pseudo[i.encode('ascii', 'ignore')] = pps[-1].split('/')[-1]
        #print("pseudo: {}".format(self.pesudo))

    def convert2primitive(self, s):
        finder = SpacegroupAnalyzer(s)
        return finder.get_primitive_standard_structure()

    def __str__(self):
        out = []
        def to_str(v):
            if isinstance(v, six.string_types):
                return "'{}'".format(v)
            return v
        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.__dict__[k1]
            out.append("&{}".format(k1.upper()) )
            sub = []
            for k2 in sorted(v1.keys()):
                sub.append("  {} = {}".format(k2, to_str(v1[k2])) )
            if k1 == "system":
                sub.append("  ibrav = 0")
                sub.append("  nat = {}".format(len(self.structure)) )
                sub.append("  ntyp = {}".format(len(self.structure.composition)) )
            sub.append("/")
            out.append(",\n".join(sub) )

        out.append("ATOMIC_SPECIES")
        for k, v in self.structure.composition.items():
            out.append("  {:2} {:.4f} {}".format(k.symbol, k.atomic_mass,
                                         self.pseudo[k.symbol]) )
        out.append("ATOMIC_POSITIONS crystal")
        for site in self.structure:
            out.append("  {:2} {:.6f} {:.6f} {:.6f}".format(site.specie.symbol, site.a,
                                                site.b, site.c) )
        out.append("K_POINTS {}".format(self.kpoints_mode) )
        if self.kpoints_mode == 'automatic':
            kpt_str = ["{}".format(i) for i in self.kpoints_grid ]
            kpt_str.extend(["{}".format(i) for i in self.kpoints_shift] )
            out.append("  {}".format(" ".join(kpt_str)) )
        elif self.kpoints_mode == 'crystal':
            #out.append("%s" % i for i in self.kpoints_grid)
            #self.kpoints_grid[-1]=self.kpoints_grid[-1].rstrip()
            #out.append("".join(self.kpoints_grid))
            out.append("  {}".format(len(self.kpoints_grid)))
            out.extend(['{:16.9F} {:14.9F} {:14.9F} {:3.1F}'.format(
                        k[0], k[1], k[2], k[3]) for k in self.kpoints_grid] )
        out.append("CELL_PARAMETERS angstrom")
        for vec in self.structure.lattice.matrix:
            out.append("  {:10.6f} {:10.6f} {:10.6f}".format(vec[0], vec[1], vec[2]))
        return "\n".join(out)


    def to_file(self, filename='in'):
        """
        Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, "w") as f:
            f.write(self.__str__()+"\n")


    @classmethod
    def from_file(cls, filename):
        sections = {}
        heading = 'none'
        headers = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 
                   'CELL_PARAMETERS', 'CONSTRAINTS', 'OCCUPATIONS', 
                   'ATOMIC_FORCES']
        to_remove = ['ibrav', 'nat', 'ntyp']

        def parse_kv(k, v):
            if 'true' in v.lower():
                return (k, True)
            elif 'false' in v.lower():
                return (k, False)
            try:
                val = float(v)
                if val.is_integer():
                    return (k, int(val))
                else:
                    return (k, val)
            except ValueError:
                if v.find('True') != -1:
                    return (k, True)
                elif v.find('False') != -1:
                    return (k, False)
                else:
                    return (k, v)
            

        def complete_subgroup(heading):
            if 'ATOMIC_SPECIES' in heading:
                global pseudo
                atom_pseudo = dcopy(sub)
                pseudo = {}
                for ppf in atom_pseudo:
                    pseudo[ppf[0]] = ppf[-1]
            if 'ATOMIC_POSITIONS' in heading:
                global species
                global coords
                atoms = dcopy(sub)
                species = [a[0] for a in atoms]
                coords = [[float(a[1]), float(a[2]), float(a[3])] for a in atoms]
            if 'K_POINTS' in heading:
                global kpoints
                global kpoints_shift
                kpt = dcopy(sub)
                #print("kpt: {}".format(kpt))
                if kpoints_mode == 'crystal':
                    kpoints = []
                    #kpoints.append('{:>5}\n'.format(kpt[0][0]))
                    kpoints.extend(['  {:11}  {:11}  {:11}   {:3}\n'.format(
                                     k[0], k[1], k[2], k[3]) for k in kpt[1:]] 
                                  )
                    kpoints_shift=None
                elif kpoints_mode == 'automatic':
                    kpoints = [int(k) for k in kpt[0][:3]]
                    kpoints_shift = [int(k) if type(k) is int else float(k) 
                            for k in kpt[0][3:]]
                #print("kps: {}".format(kpoints))
                #print("kps_shift: {}".format(kpoints_shift))
            

            if 'CELL_PARAMETERS' in heading:
                global lattice
                cell = dcopy(sub)
                lattice = [[float(v) for v in v1] for v1 in cell]
            if 'CONSTRAINTS' in heading:
                global constraints
                constraints = dcopy(sub)

        with open(filename, 'r') as f:
            lines = f.readlines()
        for line in lines:
            l = line.strip().replace("'", '').replace(',', '')
            #print("l: {}".format(l))
            if '&' in l:
                section = l[1:].lower()
                sections[section] = {}
                continue

            l = l.split()
            #print("key: {}\t\tVal: {}".format(l[0], l[-1]))
            #print("{}".format(l[-1].lower()))
            if l[0] in headers:
                if 'K_POINTS' in l[0]:
                    kpoints_mode = l[-1].strip()
                if heading != 'none':
                    complete_subgroup(heading) 
                heading = l[0]
                sub = []
                continue

            if heading != 'none':
                sub.append([w for w in l])
            elif l[0] != '/':
                k,v = parse_kv(l[0], l[-1].replace(',',''))
                sections[section][k] = v
            """elif 'true' in l[-1].lower():
                print("under true heading\n\n")
                sections[section][l[0]] = True
            elif 'false' in l[-1].lower():
                print("under false heading\n\n")
                sections[section][l[0]] = False
            elif l[0] != '/' and l[0] != 'pseudo_dir':
                print("under non pseudo heading\n\n")
                sections[section][l[0]] = (l[-1].replace("'",'').replace(',', '')) 
            elif l[0] != '/':
                print("under other heading\n\n")
                sections[section][l[0]] = l[-1][1:-2]"""
        complete_subgroup(heading)
        sections['system'] = {k:v for k,v in sections['system'].items()
                                if k not in to_remove}
        
        structure = Structure(lattice, species, coords)
        pseudo_dir=sections['control']['pseudo_dir']
        return cls(structure, pseudo_dir=pseudo_dir, control=sections['control'], 
                        system=sections['system'], electrons=sections['electrons'], 
                        ions=sections['ions'], cell=sections['cell'],
                        kpoints_mode=kpoints_mode, kpoints_grid=kpoints, 
                        kpoints_shift=kpoints_shift)


class QeMFInput(MSONable):
    '''
    Class for creating multiple input files for Quantum Espresso.  This class
    utilizes PwInput class for creating input files.
    '''
    def __init__(self, structure, pseudo_dir=None, kpoints_coarse=[1,1,1], 
                kpoints_fine=[4,4,4], kpoints_shift=[1,1,1], 
                tasks=['scf', 'wfn'], reduce_structure=False, 
                control={}, system={}, electrons={}, ions={}, 
                cell={}, config_file=None):
        """
        Initializes multiple PWSCF input files for use in WorkFlows 
        using pw.x in QuantumEspresso.

        Args:
            structure (Structure):  Input Structure
            pseudo_dir (str):       String of the directory containing pseudopotential 
                                    files.
            kpoints_coarse (list/dict):The kpoint grid for coarse grid. Can also be a 
                                    dict with different kpoint lists for each task.
                                    Default to [1,1,1].
            kpoints_fine (list/dict):The kpoint grid for fine grid. Can also be a 
                                    dict with different kpoint lists for each task.
                                    Defaults to [4,4,4].
            kpoints_shift (list):   The shift to kpoint grid used in Espresso.
                                    Options: 
                                        [0,0,0] for no shift (grid includes Gamma point)
                                        [1,1,1] for a Monkhorst Pack shift
                                    Defaults to [1,1,1].
            tasks (list):           List of PWSCF tasks to create input files.
                                    Default to ['scf', 'wfn'].
            reduce_structure (bool):If True, structure will be converted to its 
                                    reduced primitive structure.
            config_file (str):  Configuration file from which to pull user default
                                values for control, system, electrons, and ions.  Default is 
                                to look for specified file in User's HOME directory but will 
                                also accept absolute paths for alternative locations.
        """
        # Initialize attributes that are needed for other member functions.
        self.structure = structure
        self.pseudo_dir = pseudo_dir
        self.tasks = tasks

        self.kpoints_fine = kpoints_fine
        self.kpoints_coarse = kpoints_coarse
        self.kpoints_shift = kpoints_shift
        self.reduce_structure = reduce_structure

        self.control = control
        self.system = system
        self.electrons = electrons
        self.ions = ions
        self.cell = cell

        self.config_file = config_file

        # Create Input Files for list of Tasks
        for i in tasks:
            if 'fi' in i:
                if isinstance(kpoints_fine, dict):
                    kps = kpoints_fine[i]
                else:
                    kps = kpoints_fine
                if isinstance(kps[0], list):
                    kp_mode = 'crystal'
                else:
                    kp_mode = 'automatic'

                kp_shift = [1,1,1]
            else:
                if isinstance(kpoints_coarse, dict):
                    kps = kpoints_coarse[i]
                else:
                    kps = kpoints_coarse
                if isinstance(kps[0], list):
                    kp_mode = 'crystal'
                else:
                    kp_mode = 'automatic'
                kp_shift = [1, 1, 1]

                if 'co' in i:
                    kp_shift = [0, 0, 0]

            task_dicts = {}
            if 'wfn' in i:
                task_dicts['control'] = {'calculation': 'bands'}
            else:
                task_dicts['control'] = {'calculation': i}

            d = self.as_dict()
            #with open('debug.out', 'w') as fout:
            #    pprint.pprint(d, fout)
            for j in ['control', 'system', 'electrons', 'ions', 'cell']:
                subgroup = d.get(j, {}).get(i, {})
                #print("IN QEMFINPUT: j = {}\t\ttask: {}".format(j, i))
                #print("IN QEMFINPUT: subgroup: {}".format(subgroup))
                if subgroup and j in task_dicts.keys():
                    task_dicts[j].update(subgroup)
                elif subgroup:
                    task_dicts[j] = subgroup
                else:
                    task_dicts[j] = {}

            #print("IN QEMFINPUT: task_dicts: {}".format(task_dicts))
            
            setattr(self, i, PwInput(structure, pseudo_dir=pseudo_dir, 
                                control=task_dicts['control'], 
                                system=task_dicts['system'],
                                electrons=task_dicts['electrons'],
                                ions=task_dicts['ions'], cell=task_dicts['cell'],
                                kpoints_mode=kp_mode, kpoints_grid=kps, 
                                kpoints_shift=kp_shift,
                                reduce_structure=reduce_structure,
                                config_file=config_file) )

    def to_file(self):
        # Create Directory and sub directories for input files
        if not os.path.exists('./ESPRESSO'):
            os.mkdir('./ESPRESSO')
        os.chdir('./ESPRESSO')
        
        for i in self.tasks:
            # Create directory for input file
            if not os.path.exists(i):
                os.mkdir(i)
            else:
                print("Directory already exists, input files might be overwritten.")
            os.chdir(i)
            
            # Write Input File
            fin = getattr(self, i)
            fin.to_file('in')

            # Move up a Directory to prepare for next input file
            os.chdir('../')

        # Move back to original directory
        os.chdir('../')

    '''
    def as_dict(self):
        return {'@module': self.__class__.__module__,
                '@class': self.__class__.__name__,
                'structure': self.structure.as_dict(),
                'pseudo_dir': self.pseudo_dir,
                'kpoints_coarse': self.kpoints_coarse,
                'kpoints_fine': self.kpoints_fine, 
                'kpoints_shift': self.kpoints_shift,
                'tasks': self.tasks,
                'reduce_structure': self.reduce_structure,
                'control': self.control,
                'system': self.system,
                'electrons': self.electrons,
                'ions': self.ions,
                'cell': self.cell,
                'config_file': self.config_file}

    @classmethod    
    def from_dict(cls, d):
        structure = Structure.from_dict(d['structure'])
        return cls(structure, pseudo_dir=d['pseudo_dir'],
                    kpoints_coarse=d['kpoints_coarse'],
                    kpoints_fine=d['kpoints_fine'],
                    kpoints_shift=d['kpoints_shift'],
                    tasks=d['tasks'], 
                    reduce_structure=d['reduce_structure'],
                    control=d['control'], system=d['system'],
                    electrons=d['electrons'], ions=d['ions'],
                    cell=d['cell'], config_file=d['config_file'] )
    '''
                    

class PWInput(object):
    """
    Base input file class. Right now, only supports no symmetry and is
    very basic.
    """

    def __init__(self, structure, pseudo, control=None, system=None,
                 electrons=None, ions=None, cell=None, kpoints_mode="automatic",
                 kpoints_grid=(1, 1, 1),kpoints_shift=(0, 0, 0)):
        """
        Initializes a PWSCF input file.

        Args:
            structure (Structure): Input structure
            pseudo (dict): A dict of the pseudopotentials to use.
            control (dict): Control parameters. Refer to official PWSCF doc
                on supported parameters. Default to {"calculation": "scf"}
            system (dict): System parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            electrons (dict): Electron parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            ions (dict): Ions parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            cell (dict): Cell parameters. Refer to official PWSCF doc
                on supported parameters. Default to None, which means {}.
            kpoints_mode (str): Kpoints generation mode. Default to automatic.
            kpoints_grid (sequence): The kpoint grid. Default to (1, 1, 1).
            kpoints_shift (sequence): The shift for the kpoints. Defaults to
                (0, 0, 0).
        """
        self.structure = structure
        sections = {}
        sections["control"] = control or {"calculation": "scf"}
        sections["system"] = system or {}
        sections["electrons"] = electrons or {}
        sections["ions"] = ions or {}
        sections["cell"] = cell or {}
        for species in self.structure.composition.keys():
            if species.symbol not in pseudo:
                raise PWInputError("Missing %s in pseudo specification!")
        self.pseudo = pseudo
        self.sections = sections
        self.kpoints_mode = kpoints_mode
        self.kpoints_grid = kpoints_grid
        self.kpoints_shift = kpoints_shift

        '''
        # Get number of Bands needed for a good calculation
        self.bands = 0
        composition = self.structure.composition.as_dict()
        for element in self.pseudo.keys():
            with open(self.pseudo[element]) as fin:
                for line in fin:
                    if "z_valence" in line:
                        l = line.strip().split()
                        z_val = float(l[-1][:-1])
                        self.bands += ( 1.2 / 2 * z_val * composition[element] )
                        break

        #self.bands = math.ceil(1.2*sum(self.structure.atomic_numbers)/2)
        if 'nbnd' in sections['system'].keys() and sections['system']['nbnd'] > self.bands:
            pass
        else:
            print("Setting number of bands to: {}\n"
                    "Number of bands was not given or was less "
                    "than recommended number of bands".format(self.bands))
            sections['system']['nbnd'] = int(self.bands)
        '''

    def __str__(self):
        out = []
        def to_str(v):
            if isinstance(v, six.string_types):
                return "'%s'" % v
            return v
        for k1 in ["control", "system", "electrons", "ions", "cell"]:
            v1 = self.sections[k1]
            out.append("&%s" % k1.upper())
            sub = []
            for k2 in sorted(v1.keys()):
                sub.append("  %s = %s" % (k2, to_str(v1[k2])) )
            if k1 == "system":
                sub.append("  ibrav = 0")
                sub.append("  nat = %d" % len(self.structure))
                sub.append("  ntyp = %d" % len(self.structure.composition))
            sub.append("/")
            out.append(",\n".join(sub))

        out.append("ATOMIC_SPECIES")
        for k, v in self.structure.composition.items():
            out.append("  %s %.4f %s" % (k.symbol, k.atomic_mass,
                                         self.pseudo[k.symbol]))
        out.append("ATOMIC_POSITIONS crystal")
        for site in self.structure:
            out.append("  %s %.6f %.6f %.6f" % (site.specie.symbol, site.a,
                                                site.b, site.c))
        out.append("K_POINTS %s" % self.kpoints_mode)
        if self.kpoints_mode == 'automatic':
            kpt_str = ["%s" % i for i in self.kpoints_grid]
            kpt_str.extend(["%s" % i for i in self.kpoints_shift])
            out.append("  %s" % " ".join(kpt_str))
        elif self.kpoints_mode == 'crystal':
            #out.append("%s" % i for i in self.kpoints_grid)
            self.kpoints_grid[-1]=self.kpoints_grid[-1].rstrip()
            out.append("".join(self.kpoints_grid))
        out.append("CELL_PARAMETERS angstrom")
        for vec in self.structure.lattice.matrix:
            out.append("  %f %f %f" % (vec[0], vec[1], vec[2]))
        return "\n".join(out)

    def write_file(self, filename):
        """
        Write the PWSCF input file.

        Args:
            filename (str): The string filename to output to.
        """
        with open(filename, "w") as f:
            f.write(self.__str__()+"\n")

    @classmethod
    def from_file(cls, filename):
        sections = {}
        heading = 'none'
        headers = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 
                   'CELL_PARAMETERS', 'CONSTRAINTS', 'OCCUPATIONS', 
                   'ATOMIC_FORCES']
        to_remove = ['ibrav', 'nat', 'ntyp']

        def parse_kv(k, v):
            if 'true' in v.lower():
                return (k, True)
            elif 'false' in v.lower():
                return (k, False)
            try:
                val = float(v)
                if val.is_integer():
                    return (k, int(val))
                else:
                    return (k, val)
            except ValueError:
                if v.find('True') != -1:
                    return (k, True)
                elif v.find('False') != -1:
                    return (k, False)
                else:
                    return (k, v)
            

        def complete_subgroup(heading):
            if 'ATOMIC_SPECIES' in heading:
                global pseudo
                atom_pseudo = dcopy(sub)
                pseudo = {}
                for ppf in atom_pseudo:
                    pseudo[ppf[0]] = ppf[-1]
            if 'ATOMIC_POSITIONS' in heading:
                global species
                global coords
                atoms = dcopy(sub)
                species = [a[0] for a in atoms]
                coords = [[float(a[1]), float(a[2]), float(a[3])] for a in atoms]
            if 'K_POINTS' in heading:
                global kpoints
                global kpoints_shift
                kpt = dcopy(sub)
                if kpoints_mode == 'crystal':
                    kpoints = []
                    kpoints.append('{:>5}\n'.format(kpt[0][0]))
                    kpoints.extend(['  {:11}  {:11}  {:11}   {:3}\n'.format(
                                     k[0], k[1], k[2], k[3]) for k in kpt[1:]] 
                                  )
                    kpoints_shift=None
                elif kpoints_mode == 'automatic':
                    kpoints = [int(k) for k in kpt[0]]
                    kpoints_shift = [int(k) for k in kpt[1]]
            

            if 'CELL_PARAMETERS' in heading:
                global lattice
                cell = dcopy(sub)
                lattice = [[float(v) for v in v1] for v1 in cell]
            if 'CONSTRAINTS' in heading:
                global constraints
                constraints = dcopy(sub)

        with open(filename, 'r') as f:
            lines = f.readlines()
        for line in lines:
            l = line.strip().replace("'", '').replace(',', '')
            #print("L: {}".format(l))
            if '&' in l:
                section = l[1:].lower()
                sections[section] = {}
                continue

            l = l.split()
            #print("key: {}\t\tVal: {}".format(l[0], l[-1]))
            #print("{}".format(l[-1].lower()))
            if l[0] in headers:
                if 'K_POINTS' in l[0]:
                    kpoints_mode = l[-1].strip()
                if heading != 'none':
                    complete_subgroup(heading) 
                heading = l[0]
                sub = []
                continue

            if heading != 'none':
                sub.append([w for w in l])
            elif l[0] != '/':
                k,v = parse_kv(l[0], l[-1].replace(',',''))
                sections[section][k] = v
            """elif 'true' in l[-1].lower():
                print("under true heading\n\n")
                sections[section][l[0]] = True
            elif 'false' in l[-1].lower():
                print("under false heading\n\n")
                sections[section][l[0]] = False
            elif l[0] != '/' and l[0] != 'pseudo_dir':
                print("under non pseudo heading\n\n")
                sections[section][l[0]] = (l[-1].replace("'",'').replace(',', '')) 
            elif l[0] != '/':
                print("under other heading\n\n")
                sections[section][l[0]] = l[-1][1:-2]"""
        complete_subgroup(heading)
        sections['system'] = {k:v for k,v in sections['system'].items()
                                if k not in to_remove}
        
        structure = Structure(lattice, species, coords)
        return PWInput(structure, pseudo, control=sections['control'], 
                        system=sections['system'], electrons=sections['electrons'], 
                        ions=sections['ions'], cell=sections['cell'],
                        kpoints_mode=kpoints_mode, kpoints_grid=kpoints, 
                        kpoints_shift=kpoints_shift)


class PwInputError(BaseException):
    pass


class PwPpInput(abc.ABCMeta, MSONable):
    #TODO   Look into MRO and inheritance vs python version
    #       maybe use six for portability across versions like VASP?
    '''
    Docstring
    '''

    def as_dict(self):
        d = MSONable.as_dict(self)
        return d

    @abc.abstractmethod
    def write_input(self, filename):
        """Write Input Files"""
        pass

class PwBandsInput(PwPpInput):
    """
    Docstring
    """
    def __init__(self, bands=None):
        self.bands = bands


    def __str__(self):
        return

    def write_input(self, filename):
        pass

class Pw2BgwInput(MSONable):
    #TODO: Rewrite to_file() to use new self.input_pw2bgw structure instead of sections
    '''
    Base Input file Class for PW2BGW post-processing.
    '''

    def __init__(self, structure, pw2bgw_input=None, kpoints=None, 
                kpoints_shift=None, qshift=None ):
        '''
        Initializes a PW2BGW input file.

        Args:
            structure (Structure): Input Stucture.
            pw2bgw_input (dict): Input parameters.  Refer to official PW2BGW 
                Quantum Espresso doc for supported parameters. Defaults to
                {'prefix': Structure.formula, 'wfng_flag': True,
                'wfng_kgrid': False}
        '''

        self.structure = structure
        self.input_pw2bgw = pw2bgw_input or {
                                'prefix': self.structure.composition.reduced_formula,
                                'wfng_flag': True, 'wfng_kgrid': False
                                }
        self.kpoints=kpoints
        self.kpoints_shift=kpoints_shift
        self.qshift=qshift

        if self.kpoints:
            self.input_pw2bgw['wfng_kgrid'] = True
            self.input_pw2bgw['wfng_nk1'] = int(self.kpoints[0])
            self.input_pw2bgw['wfng_nk2'] = int(self.kpoints[1])
            self.input_pw2bgw['wfng_nk3'] = int(self.kpoints[2])

            self.kpoints_shift = kpoints_shift if kpoints_shift else "Monkhorst-Pack"
            self.qshift = qshift if qshift else [0, 0, 0]

            if isinstance(self.kpoints_shift, list):
                self.input_pw2bgw['wfng_dk1'] = self.kpoints_shift[0]
                self.input_pw2bgw['wfng_dk2'] = self.kpoints_shift[1]
                self.input_pw2bgw['wfng_dk3'] = self.kpoints_shift[2]
            elif "monkhorst" in self.kpoints_shift.lower(): 
                self.input_pw2bgw['wfng_dk1'] = 0.5 + self.qshift[0]*self.kpoints[0]
                self.input_pw2bgw['wfng_dk2'] = 0.5 + self.qshift[1]*self.kpoints[1]
                self.input_pw2bgw['wfng_dk3'] = 0.5 + self.qshift[2]*self.kpoints[2]
            elif 'random' in self.kpoints_shift.lower():
                self.input_pw2bgw['wfng_dk1'] = 0.47 + self.qshift[0]*self.kpoints[0]
                self.input_pw2bgw['wfng_dk2'] = 0.37 + self.qshift[1]*self.kpoints[1]
                self.input_pw2bgw['wfng_dk3'] = 0.32 + self.qshift[2]*self.kpoints[2]
            else:
                self.input_pw2bgw['wfng_dk1'] = self.qshift[0]*self.kpoints[0]
                self.input_pw2bgw['wfng_dk2'] = self.qshift[1]*self.kpoints[1]
                self.input_pw2bgw['wfng_dk3'] = self.qshift[2]*self.kpoints[2]

    def __str__(self):
        out = []
        def to_str(v):
            if isinstance(v, six.string_types):
                return "'%s'" % v
            return v

        for k1 in ['input_pw2bgw']:
            v1 = self.input_pw2bgw
            out.append("&%s" % k1.upper())
            sub = []
            for k2 in sorted(v1.keys()):
                sub.append("   %s = %s" % (k2, to_str(v1[k2])))
            sub.append("/")
            out.append(",\n".join(sub))
        return "\n".join(out)

    def to_file(self, filename):
        '''
        Write the PW2BGW input file.

        Args:
            filename (str): The string filename to output to.
        '''
        with open(filename, 'w') as fout:
            fout.write(self.__str__()+"\n")


class QeMFPw2BgwInputs(MSONable):
    #TODO: Write in functions for taking in config_file parameters from 
    #       espresso_interface_defaults.yaml ??
    #TODO: Discuss with Dr. Kedziora if this should be moved to bgw/inputs.py considering 
    #       this is a BGW program, not Espresso executable.
    '''
    Docstring
    '''
    def __init__(self, structure, pw2bgw_input={}, mf_tasks=['scf', 'wfn'], 
                config_file=None):
        self.structure = structure
        self.pw2bgw_input = pw2bgw_input
        self.mf_tasks = mf_tasks
        self.config_file = config_file
        prefix = self.structure.composition.reduced_formula

        # Set dictionaries for input files
        for i in self.mf_tasks:
            # Read in input files from kgrid.x for Kpoint parameters for Pw2Bgw.
            # Reduces user interaction and guaranteed to get the parameters 
            #   the same as was used in pw.x calculations.
            input_file = os.path.join('ESPRESSO', i, "{}.in".format(i))
            with open(input_file, 'r') as fin:
                lines = fin.readlines()
                
            kpoints = [int(j) for j in lines[0].strip().split()]
            kpoints_shift = [float(j) for j in lines[1].strip().split()]
            qshift = [float(j) for j in lines[2].strip().split()]

            # Setup parameters for input file
            task_dict = self.pw2bgw_input.get(i)

            #d = self.as_dict()
            #with open('debug_qepw2bgw.out', 'w') as fout:
            #    pprint.pprint(d, fout)
            
            setattr(self, i, Pw2BgwInput(self.structure, pw2bgw_input=task_dict,
                                        kpoints=kpoints, kpoints_shift=kpoints_shift,
                                        qshift=qshift) )
    def to_file(self):
        # Create Directory and sub directories for input files
        if not os.path.exists('./ESPRESSO'):
            os.mkdir('./ESPRESSO')
        os.chdir('./ESPRESSO')
        
        for i in self.mf_tasks:
            if not 'scf' in i.lower():
                # Change directory for writing PostProcessing input file.
                # Directory should already exist from Writing kgrid.x files.
                os.chdir(i)
                
                # Write Input File
                fin = getattr(self, i)
                fin.to_file('pp_in')

                # Move up a Directory to prepare for next input file
                os.chdir('../')

        # Move back to original directory
        os.chdir('../')

        


'''
if __name__ == "__main__":
    o = PWOutput("../../test_files/Si.pwscf.out")
    print(o.data)
    print(o.final_energy)
'''
