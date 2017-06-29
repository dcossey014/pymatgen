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


import six, math

from pymatgen import Structure
from copy import deepcopy as dcopy
from monty.re import regrep
from collections import defaultdict


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

        # Get number of Bands needed for a good calculation
        self.bands = math.ceil(1.2*sum(self.structure.atomic_numbers))
        sections['system']['nbnd'] = int(self.bands)

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




class PWInputError(BaseException):
    pass


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
