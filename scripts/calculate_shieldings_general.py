#!/usr/bin/env python

import sys
from lshiftml.models import md_analysis_light

'''
Usage: python calculate_shieldings_general.py (xyz_file)
NOTE: The passed xyz file should have its units in ANGSTROMS. Please use the bohr_2_ang.py script for conversion.
File saved as "no_eq_simulation_ase_w_cs.xyz"
'''
xyz_file= sys.argv[1]

md_analysis_light(xyz_file, 'no_eq_simulation_ase_w_cs.xyz')

print("Shielding predicitons saved as: no_eq_simulation_ase_w_cs.xyz")
