#!/usr/bin/env python
# coding: utf-8

# -----------------------------------------------
# imports
# -----------------------------------------------
import numpy as np
from ase.io import read, write
import rascal
from rascal.representations import SphericalInvariants as SOAP
from rascal.neighbourlist.structure_manager import mask_center_atoms_by_species
from rascal.utils import get_optimal_radial_basis_hypers
from rascal.models.gaptools import calculate_features

# -----------------------------------------------
# function definitions
# -----------------------------------------------

bohr2a = 0.529177

def prep_frames(frames,species):
    for i,frame in enumerate(frames):
        frame.positions *= bohr2a
        frame.cell *= bohr2a
        frame.wrap(eps=1e-13)
        mask_center_atoms_by_species(frame,[species])
    return frames

def do_fps(x, d=0):
    if d == 0 : d = len(x)
    n = len(x)
    iy = np.zeros(d, int)
    # faster evaluation of Euclidean distance
    n2 = np.sum(x**2,axis=1)
    iy[0] = 0
    dl = n2 + n2[iy[0]] - 2* np.dot(x, x[iy[0]])
    dss = []
    for i in range(1,d):
        iy[i] = np.argmax(dl)
        nd = n2 + n2[iy[i]] - 2*np.dot(x,x[iy[i]])
        dl = np.minimum(dl, nd)
        dss.append(max(dl))
    return iy,dss

# -----------------------------------------------
# specify inputs variables
# -----------------------------------------------

# compound at hand (benzene, succinic, glycine)
compound = 'glycine'
# chemical species to be considered
species = 7
# hyper parameters
[nmax, lmax, dtrans, rc, awidth, rs_rate, rs_scale, rs_exp, zeta] = [12, 9, 0.5, 14.529, 0.806, 2.587, 3.037, 3.773, 2]
 
# -----------------------------------------------
# get FPS order -- should not need to be touched
# -----------------------------------------------

# load and prepare training data
f_train = read('../Shielding_datasets/'+compound+'_ShiftML_frames_train_FPS_w_shifts.xyz', ':')
f_train = prep_frames(f_train,species)
print('training frames loaded and prepared')

# initialise features
hypers = {'soap_type' : "PowerSpectrum",
          'interaction_cutoff': rc,
          'max_radial': 12,
          'max_angular': 9,
          'gaussian_sigma_constant': awidth,
          'gaussian_sigma_type': 'Constant',
          'cutoff_function_type' : 'RadialScaling',
          'cutoff_function_parameters': dict(
              rate=rs_rate,
              scale=rs_scale,
              exponent=rs_exp
          ),
          'cutoff_smooth_width': 0.5,
          'radial_basis': 'GTO',
          'normalize' : False}
hypers = get_optimal_radial_basis_hypers(hypers, f_train, expanded_max_radial=16)
print('hypers initialised')

# initalise SOAPs
soap, manager_train = calculate_features(f_train, hypers, auto_wrap=False)
print('manager initialised')

# FPS sort the environments
fps_indices, fps_distances = do_fps(manager_train.get_features(soap))
np.savetxt('../FPS_orders/env_fps_ids_'+compound+'_'+str(species)+'_rebuilt.dat',fps_indices, fmt='%i')
