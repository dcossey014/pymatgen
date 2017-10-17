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

from fireworks import Firework, LaunchPad
from pymatgen.util.io_utils import clean_lines, micro_pyawk
from pymatgen.core.structure import Structure
from pymatgen.core.units import unitized
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.bandstructure import Spin, BandStructureSymmLine
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.bgw.kgrid import Generate_Kpath
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
        

class BgwRun(MSONable):
    def __init__(self, dir_name):
        self.dirname = dir_name
        self.output_filename = self._find_outputs(self.dirname)
        with zopen(self.output_filename, 'rt') as f:
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

        if "Sigma" in self.runtype:
            self._parse_ch_convergence() 

        self.absorption = {}
        if "Absorption" in self.runtype: 
            self.dirname = os.path.dirname(os.path.abspath(self.output_filename))
            if 'absorption_eh.dat' in os.listdir(self.dirname):
                self._parse_absorption(os.path.join(self.dirname,
                    'absorption_eh.dat'))
            if 'absorption_noeh.dat' in os.listdir(self.dirname):
                self._parse_absorption(os.path.join(self.dirname,
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

        # Parse Input of run
        inp = BgwInput.from_directory(self.dirname)
        self.inp_params = inp.as_dict()
        drop_line = self.inp_params.pop('_fw_name')
    
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
        kpt = [k.replace('.', ',') for k in kpt]
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
                key = self.key_check(m.group(1).strip())
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


    def key_check(self, key):
        return key if not '.' in key else key.replace('.', ',')


    def as_dict(self):
        d = {'BGW Version': self.ver,
                'Revision': self.rev,
                'Complex/Real': self.cmplx_real}
        d['Memory Usage'] = {'Required Memory': self.mem_req,
                            'Available Memory': self.mem_avail}
        d['Input'] = self.inp_params
        d['Output'] = {'Timings': self.timings}
        output = d['Output']
        output['Band Info'] = {'Highest Occupied Band': self.occ_band_max,
                'Max Valence Band': self.val_max_nrg,
                'Min Conduction Band': self.cond_min_nrg,
                'Fermi Energy': self.fermi_nrg,
                'Units': self.val_units}
        if "Sigma" in self.runtype:
            output['Sigma Band Avgs'] = self.band_data
            output['CH Convergence'] = self.ch_convergence

        if "Absorption" in self.runtype:
            output['Dielectric Functions'] = self.absorption
               
        return {self.runtype: d}


    def _parse_ch_convergence(self):
        i = 1
        self.ch_convergence, d = {}, {}

        with open('ch_converge.dat') as fin:
            for line in fin.readlines():
                l = line.strip().split()
                if "# k =" in line:
                    if i != 1:
                        # Update ch_convergence attribute with previous 
                        # K-point data
                        d[name]['DATA'] = data
                        self.ch_convergence.update(d)

                    # Reset Dictionary for new K-point
                    d = {}
                    name = "k-point_{}".format(i)
                    d[name] = {}
                    d[name]['K-POINT'] = ' '.join(l[-5:-2])
                    d[name]['UNITS'] = 'eV'

                    # Reset Data List for new K-point
                    data = []

                    # Increment i for new K-point
                    i += 1

                elif "nbands" in line:
                    d[name]['KEY'] = l[1:]

                elif "#" not in line:
                    data.append(l)

            # Input last k-point into ch_convergence
            d[name]['DATA'] = data
            self.ch_convergence.update(d)


    def _find_outputs(self, dirname):
        out_files = glob.glob(os.path.join(dirname, "OUT.*"))

        if len(out_files) == 1:
            return out_files[0]
        else:
            raise BgwParserError(
                    "Found {} output files for this directory: {}".format(
                        len(out_files), dirname),
                    {'err': 'Too few/many Output Files', 
                    'directory': dirname} )
        



class QeBgwRun(MSONable):
    def __init__(self, dirname='.'):
        self.dirname = dirname
        self.bgw_types = ['epsilon', 'sigma', 'kernel', 'absorption']
        esp_dirs = []
        bgw_outputs = []
        for root, dirs, files in os.walk(self.dirname):
            if "ESPRESSO" in dirs:
                esp_dirs.append(os.path.join(root, "ESPRESSO"))

            out_file = glob.glob(os.path.join(root, "OUT*"))
            if out_file:
                bgw_outputs.append(out_file[0])

        try:
            lp = Launchpad.from_file(os.path.join(os.environ['HOME'], '.fireworks', 'my_launchpad.yaml'))
            for i,esp_run in  enumerate(esp_dirs):
                #TODO: remove esp_run from esp_dirs???  Look Into this, might not be necessary
                fw_id = self._get_fw_id(os.path.dirname(esp_run))
                wf = lp.get_wf_summary_dict(fw_id)
                wf_states = wf['states']
                wf_dirs = wf['launch_dirs']
                incomplete = []
                for rtype, state in wf_states.items():
                    if state != "COMPLETED":
                        incomplete.append("{}: {}".format(rtype, wf_states[rtype]))
                if incomplete:
                    print("Incomplete runs:")
                    print('\n'.join(i for i in incomplete))
                    print("Will not compile the information at this time")
                    continue

                self.esp_data = EspressoRun(esp_run)
                #TODO: remove each BGWrun from bgw_outputs
                # Remove esp_dir from wf_dirs such that we don't scan it as 
                # bgw_run
                bgw_wf_runs = []
                for k,v in wf_dirs.items():
                    v = v[0]
                    edir = esp_run.split('/')[-2]
                    if not edir in v:
                        bgw_wf_runs.append(v)
                    else:
                        continue

                for i,bgw_run in wf_dirs.items():
                    pass

            #for fout in bgw_outputs:
            #    fw_id = self._get_fw_id(os.path.dirname(fout))
            #    linked_runs = lp.get_wf_summary_dict(fw_id)['launch_dirs'] 
        except:
            print("Unable to contact LaunchPad Server to verify "
                    "links between runs.  Continuing as if directory "
                    "contains only runs from single Molecular System.")
            if bgw_outputs:
                self.bgw_runs = self._parse_bgw_outputs(bgw_outputs)


    def _get_fw_id(self, dir):
        fw = Firework.from_file(os.path.join(dir, 'FW.json'))
        return fw.fw_id


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

