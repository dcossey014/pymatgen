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


import six, math, os, glob, re, itertools
import logging, warnings
from io import StringIO
from collections import defaultdict
import xml.etree.cElementTree as ET

from six.moves import map, zip
from six import string_types

import numpy as np

from monty.io import zopen, reverse_readfile
from monty.re import regrep
from monty.json import jsanitize, MSONable

from pymatgen.util.io_utils import clean_lines, micro_pyawk
from pymatgen.core.units import unitized
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.bandstructure import Spin, BandStructureSymmLine
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.bgw.kgrid import Generate_Kpath
from pymatgen import Structure

from copy import deepcopy as dcopy


logger = logging.getLogger(__name__)

class XmlListConfig(list):
    def __init__(self, aList):
        logger.debug("aList: {}".format(aList))
        for element in aList:
            if element:
                # treat like dict
                if len(element) == 1 or element[0].tag != element[1].tag:
                    self.append(XmlDictConfig(element))
                # treat like list
                elif element[0].tag == element[1].tag:
                    self.append(XmlListConfig(element))
            if element.text:
                logger.debug("printing text of element '{}': {}".format(
                                        element, element.text.strip()))
                text = element.text.strip()
                if text:
                    self.append(text)


class XmlDictConfig(dict):
    '''
    Example usage:
    >>> import xml.etree.cElementTree as ET

    >>> tree = ET.parse('your_file.xml')
    >>> root = tree.getroot()
    >>> xmldict = XmlDictConfig(root)

    Or, if you want to use an XML string:

    >>> root = ET.XML(xml_string)
    >>> xmldict = XmlDictConfig(root)

    And then use xmldict for what it is... a dict.
    '''
    def __init__(self, parent_element):
        logger.debug("\n\n\nparent: {}".format(parent_element))
        if parent_element.items():
            self.update(dict(parent_element.items()))
        for element in parent_element:
            element.tag = ( element.tag if '.' not in element.tag
                    else element.tag.replace('.', '_') )
            #print("element: {}".format(element))
            if element.text:
                logger.debug("text: {}\n\n".format(element.text.strip()))
            if element:
                # treat like dict - we assume that if the first two tags
                # in a series are different, then they are all different.
                if len(element) == 1 or element[0].tag != element[1].tag:
                    aDict = XmlDictConfig(element)
                # treat like list - we assume that if the first two tags
                # in a series are the same, then the rest are the same.
                else:
                    # here, we put the list in dictionary; the key is the
                    # tag name the list elements all share in common, and
                    # the value is the list itself 
                    aDict = {element[0].tag: XmlListConfig(element)} 
                # if the tag has attributes, add those to the dict
                if element.items():
                    aDict.update(dict(element.items()))
                self.update({element.tag: aDict}) 
            # This assumes that you may have an attribute in a tag
            # along with text.  
            elif element.items() and element.text:
                bDict = {element.tag: dict(element.items())} 
                cDict = bDict[element.tag]
                text = element.text.strip()
                size = int(element.attrib.get('size', 1))
                vtype = element.attrib.get('type', 'char')
                logger.debug("text: {}\ntype: {}".format(text, vtype))
                if 'columns' in element.attrib.keys():
                    if size > int(element.attrib['columns']):
                        cDict['VALUE'] = [i.split()
                                for i in text.split('\n')]
                    else:
                        cDict['VALUE'] = [self._parse_val(i, vtype)
                                for i in text.split()]
                elif size > 1:
                    cDict['VALUE'] = [self._parse_val(i, vtype)
                                for i in text.split('\n')]
                else:
                    cDict['VALUE'] = self._parse_val(text, vtype)
                self.update(bDict)
            # this assumes that if you've got an attribute in a tag,
            # you won't be having any text. 
            elif element.items():
                self.update({element.tag: dict(element.items())}) 
            # finally, if there are no child tags and no attributes, extract
            # the text
            else:
                self.update({element.tag: element.text.strip()}) 


    def _parse_val(self, val, vtype):
        logger.debug("\n\nval: {}\nvtype: {}".format(val, vtype))
        if "int" in vtype:
            return int(val)
        if "char" in vtype:
            return val
        if "real" in vtype:
            return float(val)
        if "logical" in vtype and "T" in val:
            return True
        elif "logical" in vtype and "F" in val:
            return False
        

class EspressoRun(MSONable):
    def __init__(self, run_dir, parse_charge_density=False, 
            parse_band_data=True):
        self.run_dir = run_dir
        self.parse_charge_density = parse_charge_density
        self.parse_band_data = parse_band_data

        for root, dirs, files in os.walk(run_dir, topdown=True):
            if 'scf' in dirs:
                self.espresso_runs = list(dirs)
            if 'scf' in root.split('/')[-1].lower():
                self.scf_dir = root
                self.save_dir = [s for s in dirs if '.save' in s]
                self.chg_file = "{}/{}/charge-density.xml".format(
                            root, self.save_dir[0])
                self.data_file = "{}/{}/data-file.xml".format(
                            root, self.save_dir[0]) 
                self.structure = self._parse_structure(self.scf_dir)
                break
        tree = ET.parse(self.data_file)
        xmlroot = tree.getroot()
        self.data = XmlDictConfig(xmlroot)
        
        if self.parse_charge_density:
            self.charge_dict = self._parse_chg_dens(self.run_dir)

        if self.parse_band_data:
            self.kpoints = {}
            self.band_data = self._parse_band_data(self.run_dir)
            self.band_data = self._collate_band_data(self.band_data)
            
    @property
    def efermi(self):
        return self.data['BAND_STRUCTURE_INFO']['FERMI_ENERGY']['VALUE']*27.2114

    @property
    def rec_lattice(self):
        return self.structure.lattice.reciprocal_lattice
        
    def as_dict(self):
        d = {'structure': self.structure.as_dict(),
             'data_file': self.data,
             'band_data': {'kpoints': self.kpoints,
                 'eigenvalues': self.band_data}
             }
        if self.parse_charge_density:
            d['charge_density'] = self.charge_dict['CHARGE-DENSITY']
        return d

    def gen_kpath(self, run):
        kps = self.kpoints[run]
        self.kpath = Generate_Kpath(self.structure, len(kps)-1)
        return self.kpath

    def plot_bands(self, run, filename, ylim=None, 
                    zero_to_efermi=True, usetex=False, smooth=False):
        kps = self.kpoints[run]
        latt = self.rec_lattice
        efermi = self.efermi
        evals = {Spin.up: self.band_data[run]['SORTED']}
        
        kpath = self.gen_kpath(run)
        labels = kpath.pcoords

        self.bandstructure = BandStructureSymmLine(kpoints=kps, 
                        eigenvals=evals,lattice=latt, efermi=efermi, 
                        labels_dict=labels, structure=self.structure)
        plotter = BSPlotter(self.bandstructure)
        plt = plotter.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi,
                            smooth=smooth, usetex=usetex)
        plt.savefig(filename, format=filename.split('.')[-1])


    def _parse_structure(self, d):
        with open(os.path.join(d, 'in')) as f:
            lines = f.readlines()
        for i,line in enumerate(lines):
            l = line.strip()
            if 'nat' in l:
                natoms = int(l.split()[-1].replace(',',''))
            if 'ATOMIC_POSITIONS' in l:
                apos = i + 1
                atype = l.split()[-1]
            if 'K_POINTS' in l and 'crystal' in l:
                nkps = int(lines[i+1].strip())
                kpos = i + 2
            if 'CELL_PARAMETERS' in l:
                units = l.split()[-1]
                lpos = i + 1
        #TODO : finish structure code
        #     : add KPOINTS list for use in Band structure. 
        #     : move KPOINTS to _parse_band_structure
        species = []
        atom_pos = []
        lattice = []
        #kpoints = []
        for l in lines[apos:apos+natoms]:
            l = l.strip().split()
            species.append(l[0])
            atom_pos.append([float(k) for k in l[1:4]])
        for l in lines[lpos:lpos+3]:
            l = l.strip().split()
            lattice.append([float(k) for k in l])
        #for l in lines[kpos:kpos+nkps]:
        #    l = l.strip().split()
        #    kpoints.append([float(k) for k in l[:-1]])
        s = Structure(lattice, species, atom_pos)
        return s

    def _parse_kpoints(self, rdir):
        with open(os.path.join(rdir, 'in')) as fin:
            lines = fin.readlines()
        kps = []
        for i,line in enumerate(lines):
            l = line.strip()
            if 'K_POINTS' in l and 'crystal' in l:
                nkps = int(lines[i+1].strip())
                kpos = i + 2
                break
        for l in lines[kpos:kpos+nkps]:
            l = l.strip().split()
            kps.append([float(k) for k in l[:-1]])
        return kps

    def _parse_band_data(self, run_dir):
        band_data = {}
        for root, dirs, files in os.walk(run_dir):
            if root == run_dir: 
                espr_runs = list(dirs)

            # possible solution: glob (might be faster than comparing 
            # strings across entire dir list)
            save_dir = [s for s in dirs if ".save" in s]
            if save_dir:
                run_type = root.split('/')[-1].lower()
                band_data[run_type] = {}
                band_data[run_type]['RAW'] = {}
                self.kpoints[run_type] = self._parse_kpoints(root)

            eigen_file = [s for s in files if "eigenval.xml" in s]
            if eigen_file:
                kpt = root.split('/')[-1]
                tree = ET.parse('{}/eigenval.xml'.format(root))
                xmlroot = tree.getroot()
                band_data[run_type]['RAW'][kpt] = XmlDictConfig(xmlroot)
        return band_data

    def _collate_band_data(self, bdata):
        for run_type in bdata:
            if 'SORTED' not in bdata[run_type].keys():
                num_bands = int(bdata[run_type]['RAW']['K00001']['INFO']['nbnd'])

                # Pymatgen Band Plotter needs Band Data in {'Spin.up': [[], ... , []],
                # 'Spin.down': [[], ... , []]}.  The first index (M) of the 
                # M x N array [[]] refers to the band and the second 
                # index (N) refers to the Kpoint.  If band structure is not 
                # spin polarized, only store one data set under "Spin.up".  
                # This function assumes only non-polarized band data is given.
                bdata[run_type]['SORTED'] = bds = [
                        [] for _ in range(num_bands) ]

            for kpt in sorted(bdata[run_type]['RAW']):
                kp_data = bdata[run_type]['RAW'][kpt]['EIGENVALUES']['VALUE']
                for i,v in enumerate(kp_data):
                    bds[i].append(v*27.2114)
        return(bdata)

    def _parse_chg_dens(self, run_dir):
        chg_dict = {}
        tree = ET.parse(self.chg_file)
        xmlroot = tree.getroot()
        chg_dict = XmlDictConfig(xmlroot)
        return chg_dict


class PWOutput(object):

    patterns = {
        "energies": "total energy\s+=\s+([\d\.\-]+)\sRy",
        "ecut": "kinetic\-energy cutoff\s+=\s+([\d\.\-]+)\s+Ry",
        "lattice_type": "bravais\-lattice index\s+=\s+(\d+)",
        "celldm1": "celldm\(1\)=\s+([\d\.]+)\s",
        "celldm2": "celldm\(2\)=\s+([\d\.]+)\s",
        "celldm3": "celldm\(3\)=\s+([\d\.]+)\s",
        "celldm4": "celldm\(4\)=\s+([\d\.]+)\s",
        "celldm5": "celldm\(5\)=\s+([\d\.]+)\s",
        "celldm6": "celldm\(6\)=\s+([\d\.]+)\s",
        "nkpts": "number of k points=\s+([\d]+)"
    }

    def __init__(self, filename):
        self.filename = filename
        self.data = defaultdict(list)
        self.read_pattern(PWOutput.patterns)
        for k, v in self.data.items():
            if k == "energies":
                self.data[k] = [float(i[0][0]) for i in v]
            elif k in ["lattice_type", "nkpts"]:
                self.data[k] = int(v[0][0][0])
            else:
                self.data[k] = float(v[0][0][0])

    def read_pattern(self, patterns, reverse=False,
                     terminate_on_match=False, postprocess=str):
        """
        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": "energy\(sigma->0\)\s+=\s+([\d\-\.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": "energy\(sigma->0\)\s+=\s+([\d\-\.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned
            values are lists of lists, because you can grep multiple
            items on one line.
        """
        matches = regrep(self.filename, patterns, reverse=reverse,
                         terminate_on_match=terminate_on_match,
                         postprocess=postprocess)
        self.data.update(matches)

    def get_celldm(self, i):
        return self.data["celldm%d" % i]

    @property
    def final_energy(self):
        return self.data["energies"][-1]

    @property
    def lattice_type(self):
        return self.data["lattice_type"]


if __name__ == "__main__":
    o = PWOutput("../../test_files/Si.pwscf.out")
    print(o.data)
    print(o.final_energy)
