from __future__ import division, unicode_literals, print_function

import os, json
import glob
import re
import math
import itertools
from io import StringIO
import logging
from collections import defaultdict
import pprint
import xml.etree.cElementTree as ET
import warnings
import cPickle as pickle

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
from pymatgen.io.bgw.inputs import BgwInput
from monty.json import MSONable

logger = logging.getLogger(__name__)

class XmlListConfig(list):
    def __init__(self, aList):
        #logger.debug("aList: {}".format(aList))
        for element in aList:
            if element:
                # treat like dict
                if len(element) == 1 or element[0].tag != element[1].tag:
                    self.append(XmlDictConfig(element))
                # treat like list
                elif element[0].tag == element[1].tag:
                    self.append(XmlListConfig(element))
            if element.text:
                #logger.debug("printing text of element '{}': {}".format(
                                        #element, element.text.strip()))
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
        #logger.debug("\n\n\nparent: {}".format(parent_element))
        if parent_element.items():
            self.update(dict(parent_element.items()))
        for element in parent_element:
            element.tag = ( element.tag if '.' not in element.tag
                    else element.tag.replace('.', '_') )
            #print("element: {}".format(element))
            if element.text:
                #logger.debug("text: {}\n\n".format(element.text.strip()))
                pass
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
                #logger.debug("text: {}\ntype: {}".format(text, vtype))
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
        #logger.debug("\n\nval: {}\nvtype: {}".format(val, vtype))
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
            if line.find("BerkeleyGW branch") == 0:
                logger.debug("parsing version: {}".format(line))
                self._parse_version(line)
            if line.find("version, run") != -1:
                logger.debug("parsing runtype: {}".format(line))
                self._parse_runtype(line) 
            if "MB per PE" in line:
                self._parse_memory(line) 
            if "grid)" in line:
                self._parse_band_info(line) 
            if "CPU (s)" in line:
                self._parse_timings(i, lines) 
            if line.find("number of bands") == 0:
                self.num_bands = int(line.split()[-1])

            if "Sigma" in self.runtype and "Number of bands" in line:
                    self.sigma_bands_calulated = int(line.split()[-1])

            if "Sigma" in self.runtype and "Symmetrized values" in line:
                    self._parse_sigma_band_avg(i, self.sigma_bands_calulated, lines) 

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
            #logger.debug('setting available memory: mem_avail: {}'.format(
            #                self.mem_avail))
        else:
            self.mem_req += float(l.split()[-4])
            #logger.debug('setting memory required: {}'.format(self.mem_req))

    def _parse_runtype(self, stream):
        l = stream.split()
        print("In _parse_runtype: {}".format(l))
        logger.debug("In _parse_runtype: {}".format(l))
        self.runtype = l[0]
        self.cmplx_real = l[2]
        logger.debug("complx/real: {}".format(self.cmplx_real))
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
        logger.debug("in Sig Band AVg: i: {}\tj: {}".format(i,j))
        logger.debug("parsing kpt: {}".format(stream[i+2]))
        kpt = stream[i+2].strip().split()[2:5]
        kpt = [k.replace('.', ',') for k in kpt]
        kpt_str = '   '.join(kpt)
        logger.debug("kpt string: {}".format(kpt_str))
        key = stream[i+4].strip().split()
        logger.debug("Key: {}".format(key))
        d = {}
        for line in stream[i+5 : i+j+5 ]:
            l = line.strip().split()
            logger.debug("parsing values: {}".format(l))
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

    def to_file(self, filename):
        with open(filename, 'w') as fout:
            pprint.pprint(self.as_dict(), fout, indent=1 )

    def _parse_ch_convergence(self):
        i = 1
        self.ch_convergence, d = {}, {}

        with open(os.path.join(self.dirname, 'ch_converge.dat')) as fin:
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
                    d[name]['K-POINT'] = ' '.join(l[-4:-3])
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

