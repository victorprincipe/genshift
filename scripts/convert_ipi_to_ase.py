#!/bin/python
import numpy as np                                                                                                                                                                                  
from ase.io import read,write                                                                                                                                                                       
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="input i-PI style extended xyz", type=str)
parser.add_argument("outfile", help="output ASE style extended xyz", type=str)
args = parser.parse_args()

infile = args.infile
outfile = args.outfile
traj = read(infile,':')                                                                                                                                                              
cells = [] 
with open(infile) as f: 
    for il,l in enumerate(f): 
        if '#' in l: 
            c = l.strip().split(' ') 
            c = np.array([x for x in c if x.strip()][2:8]) 
            c = np.array([float(x) for x in c]) 
            cells.append(c) 
for ifrm,frm in enumerate(traj): 
    frm.set_pbc(True) 
    cell = cells[ifrm] 
    frm.set_cell(cell) 
    pos = frm.get_scaled_positions() 
write(outfile,traj)
