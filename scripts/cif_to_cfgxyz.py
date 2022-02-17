'''
Convert a .cif file to .xyz and .cfg. Also removes unnecessary array from the xyz and cfg files.
'''

from sys import argv
from ase.io import read,write

file = argv[1]

def cif_to_cfgxyz(cif_file):
    cif = read(cif_file)
    cif.set_array('spacegroup_kinds', None)
    write(cif_file[:-3]+'xyz', cif)
    write(cif_file[:-3]+'cfg', cif)

    return

cif_to_cfgxyz(file)
