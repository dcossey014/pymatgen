from __future__ import division, unicode_literals, print_function

import os
import glob
import re
import math
import itertools
from io import StringIO
import logging
from collections import defaultdict
import xml.etree.cElementTree as ET
import warnings

from six.moves import map, zip
from six import string_types

import numpy as np

from monty.io import zopen, reverse_readfile
from monty.re import regrep
from monty.json import jsanitize

from pymatgen.util.io_utils import clean_lines, micro_pyawk
from pymatgen.core.structure import Structure
from pymatgen.core.units import unitized
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from monty.json import MSONable

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
            if 'scf' in root.split('/')[-1]:
                self.scf_dir = root
                self.save_dir = [s for s in dirs if '.save' in s]
                self.chg_file = "{}/{}/charge-density.xml".format(
                            root, self.save_dir[0])
                self.data_file = "{}/{}/data-file.xml".format(
                            root, self.save_dir[0]) 
                break
        tree = ET.parse(self.data_file)
        xmlroot = tree.getroot()
        self.data = XmlDictConfig(xmlroot)
        
        if self.parse_charge_density:
            self.charge_dict = self._parse_chg_dens(self.run_dir)

        if self.parse_band_data:
            self.band_data = self._parse_band_data(self.run_dir)

    def as_dict(self):
        d = {'data_file': self.data,
             'band_data': self.band_data}
        if self.parse_charge_density:
            d['charge_density'] = self.charge_dict['CHARGE-DENSITY']
        return d

    def _parse_band_data(self, run_dir):
        band_data = {}
        for root, dirs, files in os.walk(run_dir):
            if root == run_dir: 
                espr_runs = list(dirs)

            # possible solution: glob (might be faster than comparing 
            # strings across entire dir list)
            save_dir = [s for s in dirs if ".save" in s]
            if save_dir:
                run_type = root.split('/')[-1]
                band_data[run_type] = {}

            eigen_file = [s for s in files if "eigenval.xml" in s]
            if eigen_file:
                kpt = root.split('/')[-1]
                tree = ET.parse('{}/eigenval.xml'.format(root))
                xmlroot = tree.getroot()
                band_data[run_type][kpt] = XmlDictConfig(xmlroot)
        return band_data

    def _parse_chg_dens(self, run_dir):
        chg_dict = {}
        tree = ET.parse(self.chg_file)
        xmlroot = tree.getroot()
        chg_dict = XmlDictConfig(xmlroot)
        return chg_dict


class BgwRun(MSONable):
    def __init__(self, filename):
        self.filename = filename
        out_files = []
        with zopen(filename, 'rt') as f:
            self._parse(f)


    def _parse(self, stream):
        lines = stream.readlines()
        self.mem_req = 0.0
        self.runtype = 'unset'
        self.band_data = {}
        self.timings = {}
        for i, line in enumerate(lines):
            line = line.strip()
            if line.find("BerkeleyGW") == 0:
                self._parse_version(line)
            if line.find("version  Run") != -1:
                self._parse_runtype(line) 
            if "MB per PE" in line:
                self._parse_memory(line) 
            if "grid)" in line:
                self._parse_band_info(line) 
            if "CPU [s]" in line:
                self._parse_timings(i, lines) 
            if line.find("ndiag") == 0:
                self.num_bands = int(line.split()[-1])

            if "Sigma" in self.runtype and "Symmetrized values" in line:
                self._parse_sigma_band_avg(i, self.num_bands, lines) 

            self.absorption = {}
            if "Absorption" in self.runtype: 
                abs_dir = os.path.dirname(os.path.abspath(self.filename))
                if 'absorption_eh.dat' in os.listdir(abs_dir):
                    self._parse_absorption(os.path.join(abs_dir,
                        'absorption_eh.dat'))
                if 'absorption_noeh.dat' in os.listdir(abs_dir):
                    self._parse_absorption(os.path.join(abs_dir,
                        'absorption_noeh.dat'))

    def _parse_memory(self, stream):
        l = stream.strip()
        if l.find("Memory available:") != -1: 
            self.mem_avail = float(l.split()[-4])
            logger.debug('setting available memory: mem_avail: {}'.format(
                            self.mem_avail))
        else:
            self.mem_req += float(l.split()[-4])
            logger.debug('setting memory required: {}'.format(self.mem_req))

    def _parse_runtype(self, stream):
        l = stream.split()
        self.runtype = l[0]
        self.cmplx_real = l[1]
        if "Sigma" in self.runtype:
            self.band_avgs = {} 
    
    def _parse_version(self, stream):
        l = stream.split()
        self.ver, self.rev = (l[2], l[-1]) if len(l) < 6 else (
                                "{} {}".format(l[2], l[3]), l[-1])
        #return {"Version": ver, "Revision": rev}

    def _parse_band_info(self, stream):
        l = stream.split()
        if "Highest occupied" in stream:
            self.occ_band_max = l[-1]
        if 'Valence' in stream:
            self.val_max_nrg, self.val_units =  (l[-2], l[-1])
        if 'Conduction' in stream:
            self.cond_min_nrg, self.cond_units = (l[-2], l[-1])
        if 'Fermi' in stream:
            self.fermi_nrg, self.fermi_units = (l[-2], l[-1])

    def _parse_sigma_band_avg(self, i, j, stream):
        kpt = stream[i+2].strip().split()[2:5]
        kpt_str = '   '.join(kpt)
        key = stream[i+4].strip().split()
        d = {}
        for line in stream[i+5 : i+j+5 ]:
            l = line.strip().split()
            d[l[0]] = {val: l[k] for k,val in enumerate(key)
                    if k != 0}
        self.band_data[kpt_str] = d
            
    def _parse_timings(self, i, stream):
        for line in stream[i+2:]:
            l = line.strip()
            if self.runtype == 'Epsilon' and len(l.split()) > 3 or \
                    self.runtype == 'Sigma' and len(l.split()) > 3:
                m = re.match(
                    "([a-zA-Z ()-/_.0-9]+)[ :]+ (\d+.\d+)\s+"+\
                            "(\d+.\d+)\s+(\d+)", l)
            elif len(l.split()) > 4:
                m = re.match(
                    "([a-zA-Z ()-/_.0-9]+)[ :]+ [a-zA-Z .0()]+\s+"+\
                            "(\d+.\d+)\s+(\d+.\d+)\s+(\d+)", l)
            elif len(l.split()) == 4 and "TOTAL" in l:
                m = re.match(
                    "([a-zA-Z]+)[ :]+ \([a-zA-Z .0]+\)\s+"+\
                            "(\d+.\d+)\s+(\d+.\d+)", l)
            elif len(l.split()) == 3 and "TOTAL" in l:
                m = re.match(
                    "([a-zA-Z]+)[ :]+ (\d+.\d+)\s+"+\
                            "(\d+.\d+)", l)

            if m.groups():
                m.groups()
                key = m.group(1).strip()
                cpu_time = float(m.group(2).strip())
                wall_time = float(m.group(3).strip())
                if wall_time > 0:
                    self.timings[key] = {'CPU [s]': cpu_time,
                                'Walltime [s]': wall_time}
                    if len(m.groups()) > 3 and int(m.group(4).strip()) > 0:
                        calls = int(m.group(4).strip())
                        self.timings[key]['# of Calls'] = calls

    def _parse_absorption(self, filename):
        ftype = 'No Exciton Hole' if ( 'noeh' in 
                filename.split('.')[0].split('_')[-1] ) else (
                'Exciton Hole' )
        with open(filename, 'r') as f:
            ab_lines = f.readlines()
        
        d = {ftype: {}}
        key = []
        coords = []
        for line in ab_lines:
            l = line.strip().split()
            if "#" in line and not "Column 1" in line:
                d[ftype][l[-1]] = []
                key.append(l[-1])
            elif not "Column 1" in line and l:
                for i,k in enumerate(key, start=1):
                    d[ftype][k].append([l[0], l[i]])
        self.absorption.update(d)

    def as_dict(self):
        d = {'BGW Version': self.ver,
                'Revision': self.rev,
                'Complex/Real': self.cmplx_real}
        d['Memory Usage'] = {'Required Memory': self.mem_req,
                            'Available Memory': self.mem_avail}
        d['Runtimes'] = self.timings
        d['Band Info'] = {'Highest Occupied Band': self.occ_band_max,
                'Max Valence Band': self.val_max_nrg,
                'Min Conduction Band': self.cond_min_nrg,
                'Fermi Energy': self.fermi_nrg,
                'Units': self.val_units}
        if "Sigma" in self.runtype:
            d['Sigma Band Avgs'] = self.band_data

        if "Absorption" in self.runtype:
            d['Dielectric Functions'] = self.absorption
               
        return {self.runtype: d}


class QeBgwRun(MSONable):
    def __init__(self, dirname='.'):
        self.dirname = dirname
        self.bgw_types = ['epsilon', 'sigma', 'kernel', 'absorption']
        for root, dirs, files in os.walk(self.dirname):
            if "ESPRESSO" in dirs:
                esp_dir = os.path.join(root, "ESPRESSO")
                break
        if esp_dir:
            self.esp_runs = EspressoRun(esp_dir)
        
        bgw_outputs = []
        for root, dirs, files in os.walk(self.dirname):
            out_file = glob.glob(os.path.join(root, "OUT*"))
            if out_file:
                bgw_outputs.append(out_file[0])
        if bgw_outputs:
            self.bgw_runs = self._parse_bgw_outputs(bgw_outputs)

    def _parse_bgw_outputs(self, bgw_files):
        print("working with these outputs: {}".format(bgw_files))
        d = {}
        for i in bgw_files:
            print("working on file: {}".format(i))
            tmp_out = BgwRun(i)
            if tmp_out.runtype not in d.keys() and tmp_out.timings:
                d.update(tmp_out.as_dict())
            elif tmp_out.runtype in d.keys():
                raise BgwParserError(
                        "Found more than one output file for "
                        "Runtype: {}".format(tmp_out.runtype), 
                        {'err': 'Duplicate Run', 'file': i} )
            else:
                raise BgwParserError(                                                           
                        "Could not parse timings.  Make sure the run "
                        "completed successfully.",
                        {'err': 'Incomplete Run'} )
        return d

    def as_dict(self):
        d = {'ESPRESSO': self.esp_runs.as_dict(), 'BGW': self.bgw_runs}
        return d


class BgwParserError(Exception):
    def __init__(self, msg, err):
        emsg = []
        emsg.append("\n" + ("*" * (len(msg) + 6)))
        emsg.append("** " + msg + " **")
        emsg.append("*" * (len(msg) + 6))
        emsg = '\n'.join(i for i in emsg)
        super(BgwParserError, self).__init__(emsg)
