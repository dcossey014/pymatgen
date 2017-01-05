from kgrid import Kgrid
from pymatgen import Structure

s = Structure.from_file('/workspace/dec014/BGW/mp-33088_Cr2FeO4.cif')
tkg = Kgrid(s, kpoints=[5,5,5], offset_type='Gamma')
tkg.write_input('tmp.out')

