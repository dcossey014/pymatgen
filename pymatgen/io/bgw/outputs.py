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
from sklearn.metrics import mean_squared_error

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
    # Parses a BGW subtask output. Called from BgwDB.run_task() in 
    # custodian_job.py
    def __init__(self, runtype, dir_name):
        self.dirname = dir_name
        self.runtype = runtype
        print("Processing",runtype,"in BgwRun")
        print("gk: in BgwRun.__init__, self.dirname=",self.dirname)
        # get OUT.* file
        self.output_filename = self._find_outputs(self.dirname)
        with zopen(self.output_filename, 'rt') as f:
            self._parse(f)


    def _parse(self, stream):
        # This is where the parsing of various BGW programs is done
        lines = stream.readlines()
        self.mem_req = 0.0
        self.band_data = {}
        self.timings = {}

        def _parse_out_file(lines):
            # Parse common features of OUT.xxx

            for i, line in enumerate(lines):
                line = line.strip()
                if line.find("BerkeleyGW branch") == 0:
                    self._parse_version(line)
                if line.find("version, run") != -1: 
                    self._parse_cmplx_real(line) 
                if "MB per PE" in line:
                    self._parse_memory(line) 
                if "grid)" in line:
                    self._parse_band_info(line) 
                if "CPU (s)" in line:
                    # this has some runtype dependancies
                    self._parse_timings(i, lines) 
                if line.find("Number of bands") != -1:
                    self.num_bands = int(line.split()[-1])

        # call _parse_out_file for all, i.e. epsilon, sigma, kernel, absorption
        _parse_out_file(lines)

        if "epsilon" in self.runtype:
            print("Processing additional Epsilon data")
            # parse epsilon.log in the epsilon directory
            self._parse_epsilon_log()
            # parse ch_converged.dat in the epsilon directory
            self._parse_epsilon_conv()

        if "sigma" in self.runtype:
            self.band_avgs = {} 
            for i, line in enumerate(lines):
                line = line.strip()
                if "Symmetrized values" in line:
                    self._parse_sigma_band_avg(i, self.num_bands, lines) 
            # parse sigma ch_convergence.dat
            self._parse_ch_convergence() 

        if "kernel" in self.runtype:
            pass

        self.absorption = {}
        if "absorption" in self.runtype: 
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

    def _parse_cmplx_real(self, stream):
        l = stream.split()
        self.cmplx_real = l[2]
        if "sigma" in self.runtype:
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

    def _parse_epsilon_log(self):

        self.epsilon_log={}
        eps_fh=open(os.path.join(self.dirname,'epsilon.log'))
        
        iqpt=0
        qpts=[]
        num_g_min=99999999
        epsinv_min=1.0
        epsinv_off_max=0.0
    
        line=eps_fh.readline() 
        while line:
            while line:
                if "q=" in line:
                    ll=line.split()
                    qpt=[float(ll[1]),float(ll[2]),float(ll[3])]
                    qpts.append(qpt)
                    break
                line=eps_fh.readline() 
    
            while line:
                if 'inverse epsilon' in line:
                    break
                line=eps_fh.readline() 
    
            line=eps_fh.readline() 
            G_1=[]
            G_2=[]
            epsinv=[]
            num_g=0
            while line:
                if line in [' \n', ' \r\n', '\n', '\r\n']:
                    break
                ll=line.split()
                g1=[int(ll[0]), int(ll[1]), int(ll[2])]
                g2=[int(ll[3]), int(ll[4]), int(ll[5])]
                epsinv_r=float(ll[6])
                if len(ll) == 8:
                    epsinv_c=float(ll[7])
                else:
                    epsinv_c=0.0
                epsinv_len=math.sqrt(epsinv_r**2+epsinv_c**2)
            
                G_1.append(g1)
                G_2.append(g2)
                epsinv.append(epsinv_len)
                
                num_g=num_g+1
                line=eps_fh.readline() 
    
            line=eps_fh.readline() 
            num_g_min=min(num_g,num_g_min)
            epsinv_min=min(epsinv_len,epsinv_min)
            epsinv_off_max=max(epsinv[num_g-2],epsinv_off_max)
    
        self.epsilon_log["minimum number of G vectors"]=num_g_min
        self.epsilon_log["minimum EpsilonInverse_GmxGmx"]=epsinv_min
        self.epsilon_log["maximum EpsilonInverse_Gmx(Gmx-1)"]=epsinv_off_max
    
        eps_fh.close()
        

    def _parse_epsilon_conv(self):
    
        self.epsilon_chi_convergence={}
    
        # arrays for sorted plot
        number_bands=[]
        labels_00=[]
        labels_MM=[]
    
        #for file_name in files_list:
        chi_fh=open(os.path.join(self.dirname,'chi_converge.dat'))
        
        iqpt=0
        qpts=[]
        chi00mx=[]       #  chi_00(q,0) with maximum number of bands, G=0  
        chi00mx_extr=[]  #  chi_00(q,0) extrapolated valule with maximum number of bands  
        chi00mx_err=[]   #  chi_00(q,0) error with respect to extrapolated value
        chiMMmx=[]       #  chi_GmGm(q,0) with maximum number of bands Gm=G_max
        chiMMmx_extr=[]  #  chi_GmGm(q,0) with maximum number of bands extrapolated value
        chiMMmx_err=[]   #  chi_GmGm(q,0) with maximum number of bands error wrt extrapolated
    
        outstrs=[]
    
        line=chi_fh.readline() 
        while line:
            while line:
                if "q=" in line:
                    #print line
                    (s1,s2,qx,qy,qz,s3,s4)=line.split()
                    qpt=[float(qx),float(qy),float(qz)]
                    qpts.append(qpt)
                    break
            line=chi_fh.readline()
            line=chi_fh.readline()
            ibnd=0
            while line:
                #print line
                if line in [' \n', ' \r\n', '\n', '\r\n']:
                    #print "blank line:",line
                    break
                (i, chi00, extrp00, chiGmGm, extrpGmGm) = line.split()
                ibnd=int(i)
                line=chi_fh.readline()
            # at max number of bands. collect data for each kpoint
            chi00mx.append(float(chi00))
            chi00mx_extr.append(float(extrp00))
            chi00mx_err.append(-float(extrp00)+float(chi00))
            chiMMmx.append(float(chiGmGm))
            chiMMmx_extr.append(float(extrpGmGm))
            chiMMmx_err.append(-float(extrpGmGm)+float(chiGmGm))
            iqpt=iqpt+1
            line=chi_fh.readline()
            num_qpts=iqpt
            num_bands=ibnd
        number_bands.append(num_bands)
        outstrs.append("# number of bands = "+str(num_bands)+" used in chi\n")
        outstrs.append("# number of q-points="+str(num_qpts)+"\n")
        outstrs.append("# chi00mx,  chi00mx_extr, chi00mx_err,  chiMMmx,  chiMMmx_extr, chiMMmx_err\n")
        for iqpt in range(0,num_qpts):
            outstrs.append("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n" % (chi00mx[iqpt], chi00mx_extr[iqpt], chi00mx_err[iqpt], 
                chiMMmx[iqpt], chiMMmx_extr[iqpt], chiMMmx_err[iqpt]))
        
        chi00_rmse=math.sqrt(mean_squared_error(chi00mx_extr,chi00mx))
        chiMM_rmse=math.sqrt(mean_squared_error(chiMMmx_extr,chiMMmx))
        chi00_abse=np.array(chi00mx_err)
        chi00_maxe=np.amax(chi00_abse)
        chiMM_abse=np.array(chiMMmx_err)
        chiMM_maxe=np.amax(chiMM_abse)
    
        self.epsilon_chi_convergence["errors vs qpoints table"]=outstrs
        self.epsilon_chi_convergence["chi_00(q,0) RMSE"]=chi00_rmse
        self.epsilon_chi_convergence["chi_00(q,0) Max Error"]=chi00_maxe
        self.epsilon_chi_convergence["chi_GmGm(q,0) RMSE"]=chiMM_rmse
        self.epsilon_chi_convergence["chi_00(q,0) Max Error"]=chiMM_maxe

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
            if self.runtype == 'epsilon' and len(l.split()) > 3 or \
                    self.runtype == 'sigma' and len(l.split()) > 3:
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
        # this is where the data gets put into an overall dictionary
        # for each BGW task when called by BgwDB.get_dict()
        d = {'BGW Version': self.ver,
                'Revision': self.rev,
                'Complex/Real': self.cmplx_real}
        d['Memory Usage'] = {'Required Memory': self.mem_req,
                            'Available Memory': self.mem_avail,
                            'UNITS': 'MB per PE'}
        d['Input'] = self.inp_params
        d['Output'] = {'Timings': self.timings}
        output = d['Output']
        output['Band Info'] = {'Highest Occupied Band': self.occ_band_max,
                'Max Valence Band': self.val_max_nrg,
                'Min Conduction Band': self.cond_min_nrg,
                'Fermi Energy': self.fermi_nrg,
                'Units': self.val_units}
        if "epsilon" in self.runtype:
            output["Chi Convergence"] = self.epsilon_chi_convergence
            output["Epsilon Log"] = self.epsilon_log
        if "sigma" in self.runtype:
            output['Sigma Band Avgs'] = self.band_data
            output['Chi Convergence'] = self.ch_convergence

        if "absorption" in self.runtype:
            output['Dielectric Functions'] = self.absorption
               
        #return {self.runtype: d}
        return d


    def _parse_ch_convergence(self):
        i = 1
        self.ch_convergence, d = {}, {}

        with open(os.path.join(self.dirname, 'ch_converge.dat') ) as fin:
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
        

class BgwParserError(Exception):
    def __init__(self, msg, err):
        emsg = []
        emsg.append("\n" + ("*" * (len(msg) + 6)))
        emsg.append("** " + msg + " **")
        emsg.append("*" * (len(msg) + 6))
        emsg = '\n'.join(i for i in emsg)
        super(BgwParserError, self).__init__(emsg)

