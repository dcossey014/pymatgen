## for Surface Energy Calculation
from __future__ import division, unicode_literals

__author__ = "Zihan XU"
__version__ = "0.1"
__email__ = "vivid0036@gmail.com"
__date__ = "6/2/15"


from pymatgen.io.vasp.sets import MPVaspInputSet, DictVaspInputSet
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure


"""
for VASP Input
"""
class MPSlabVaspInputSet(DictVaspInputSet):
    """
    Class for writing a slab vasp run with LDA method.

    Args:
        user_incar_settings(dict): A dict specifying additional incar
            settings, default to None
            (ediff_per_atom=False)
        kpoints0: specify kpts[0] = kpoints0, default to []
        k_product: kpts[0][0]*a. Decide k density without kpoint0,
        default to 45
        potcar_functional: default to PBE
        bulk: default to False


        **kwargs:
            Other kwargs supported by :class:`DictVaspInputSet`.
    """
    def __init__(self, user_incar_settings=None, kpoints0=[],
                  k_product=50, potcar_functional='PBE', bulk=False, **kwargs):
        vis = MPVaspInputSet(ediff_per_atom=False).as_dict()
        DictVaspInputSet.__init__(self, "MaterialsProject Slab", vis["config_dict"],
                                  **kwargs)
        incar_settings_basic = {"NPAR": 8, "EDIFF": 1e-6, "EDIFFG": -0.01,
                                "ENCUT": 400, "ISMEAR": 0, "SIGMA": 0.05, "ISIF": 3}

        if bulk:
             self.incar_settings.update(incar_settings_basic)
        else:
            incar_settings_basic["ISIF"] = 2
            incar_settings_basic["AMIN"] = 0.01
            incar_settings_basic["AMIX"] = 0.2
            incar_settings_basic["BMIX"] = 0.001
            incar_settings_basic["NELMIN"] = 8
            self.incar_settings.update(incar_settings_basic)
        self.user_incar_settings = user_incar_settings or {}
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

        self.k_product = k_product
        self.kpoints0 = kpoints0
        self.potcar_functional = potcar_functional
        self.bulk = bulk


    def get_kpoints(self, structure):
        """
        kpoint0 is the first consideration,
        k_product is a second choice, default to 40
        """

        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter

        structure = Structure.from_dict(structure.as_dict())

        kpt = super(MPSlabVaspInputSet, self).get_kpoints(structure)
        kpt.comment = "Automatic mesh"
        kpt.style = 'Gamma'

        # use k_product to calculate kpoints, k_product = kpts[0][0] * a
        abc = structure.lattice.abc
        kpt_calc = [int(self.k_product/abc[0]+0.5),
                    int(self.k_product/abc[1]+0.5), 1]
        self.kpt_calc = kpt_calc
        # calculate kpts (c direction) for bulk. (for slab, set to 1)
        if self.bulk:
            kpt_calc[2] = int(self.k_product/abc[2]+0.5)

        # kpoint0 is prior to k_product
        if self.kpoints0:
            kpt.kpts[0] = self.kpoints0
            print "kpoint0", kpt.kpts
        else:
            kpt.kpts[0] = kpt_calc
            print 'kpt_calc:\n', kpt_calc

        return kpt

    def get_incar(self, structure):

        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter

        structure = Structure.from_dict(structure.as_dict())

        abc = structure.lattice.abc
        kpt_calc = [int(self.k_product/abc[0]+0.5),
                    int(self.k_product/abc[1]+0.5),
                    int(self.k_product/abc[1]+0.5)]

        if self.kpoints0:
            kpts = self.kpoints0
        else:
            kpts = kpt_calc

        if kpts[0]<5 and kpts[1]<5:
            if not self.bulk:
                self.incar_settings.update(
                    {"ISMEAR": 0})
            else:
                if kpts[2]<5:
                    self.incar_settings.update(
                        {"ISMEAR": 0})
        if self.user_incar_settings:
                self.incar_settings.update(self.user_incar_settings)

        incr = super(MPSlabVaspInputSet, self).get_incar(structure)

        return incr

    def as_dict(self):
        d = super(MPSlabVaspInputSet, self).as_dict()
        d.update({
            "kpoints0": self.kpoints0,
            "potcar_functional": self.potcar_functional,
            "user_incar_settings": self.user_incar_settings
        })
        return d

    def from_dict(cls, d):
        return cls(kpoints0=d.get("kpoints0", []),
                   user_incar_settings=d.get("user_incar_settings", None),
                   potcar_functional=d.get("potcar_functional", None))

    def get_all_vasp_input(self, structure):
        """
        Returns all input files as a dict of {filename: vaspio object}

        Args:
            structure (Structure/IStructure): Structure to generate vasp
                input for.
        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        """

        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter

        structure = Structure.from_dict(structure.as_dict())

        data = {'INCAR': self.get_incar(structure),
                'KPOINTS': self.get_kpoints(structure),
                'POSCAR': self.get_poscar(structure),
                'POTCAR': self.get_potcar(structure)}
        return data


def get_input_mp(element, miller_index, api_key, min_slab_size=10, min_vacuum_size=10,
                 symprec=0.001, angle_tolerance=5):
    """
    element: str, element name of Metal
    miller_index: hkl, e.g. [1, 1, 0]
    api_key: to get access to MP DB
    """
    # This initializes the REST adaptor. Put your own API key in.
    mprest = MPRester(api_key)
    #first is the lowest energy one
    prim_unit_cell = mprest.get_structures(element)[0]
    spa = SpacegroupAnalyzer(prim_unit_cell,  symprec=symprec,
                             angle_tolerance=angle_tolerance)
    conv_unit_cell = spa.get_conventional_standard_structure() 
    slab_metal = SlabGenerator(conv_unit_cell, miller_index, min_slab_size, min_vacuum_size, 
                               lll_reduce=False, center_slab=False, primitive=False)
    oriented_u_cell = slab_metal.oriented_unit_cell
    slab_cell = slab_metal.get_slab()

    return oriented_u_cell, slab_cell

