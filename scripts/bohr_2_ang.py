import sys
import ase
from ase.io import read, write

'''
Usage: python bohr_2_ang.py (xyz file)
Convert an ASE-formatted xyz trajectory from bohr to angstroms. This is required for Matthias's general shielding model.
New file is saved as (old_filename)_ang.xyz
'''

def bohr_2_ang(file):
    bohr = read(file, index=":") # Read file

    for frame in bohr:
        frame.cell *= ase.units.Bohr  # Convert cell params to Bohr
        frame.positions *= ase.units.Bohr  # Convert positions to Bohr

    write(file[:-4]+'_ang.xyz', bohr) # Save file 
    print("Conversion complete. File saved as: "+file[:-4]+'_ang.xyz' )

bohr_2_ang(sys.argv[1])

