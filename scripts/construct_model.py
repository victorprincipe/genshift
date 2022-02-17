#!/usr/bin/env python
# coding: utf-8

# -----------------------------------------------
# imports
# -----------------------------------------------
import numpy as np
from ase.io import read, write
import pickle
import rascal
from rascal.representations import SphericalInvariants as SOAP
from rascal.neighbourlist.structure_manager import mask_center_atoms_by_species
from rascal.utils import get_optimal_radial_basis_hypers
from rascal.models import KRR
from rascal.models.sparse_points import SparsePoints
from rascal.models.gaptools import calculate_features, build_sparse_list, compute_kernels, fit_gap_simple

# -----------------------------------------------
# specify inputs variables
# -----------------------------------------------
# compound at hand (benzene, succinic, glycine)
compound = 'glycine'
# chemical species to be considered
species  = 7
# hyper parameters
[nmax, lmax, dtrans, rc, awidth, rs_rate, rs_scale, rs_exp, zeta] = [12, 9, 0.5, 7.461, 0.333, 1.860, 2.159, 5.590, 2]
# number of subsampling models
nmodels  = 16
# fraction of (FPS ordered) environments to include in the active set
fraction = 0.25
 
# -----------------------------------------------
# function definitions -- should not need to be touched
# -----------------------------------------------
def prep_frames(frames,species):
    bohr2ang = 0.529177
    chemical_shifts = []
    number_of_atoms = np.zeros(len(frames), int)
    for i,frame in enumerate(frames):
        # convert from atomic units to Angstrom
        frame.positions *= bohr2ang
        frame.cell *= bohr2ang
        # wrap all atoms into unit cell
        frame.wrap(eps=1e-13)
        # mask all atoms except chemical species of interest to only calculate atomic descriptors for the latter
        mask_center_atoms_by_species(frame,[species])
        mask_species = np.where(frame.get_atomic_numbers() == species)[0]
        number_of_atoms[i] = len(mask_species)    
        chemical_shifts += list(frame.arrays['CS_total'][mask_species])
    return frames, np.asarray(chemical_shifts), number_of_atoms

def get_rmse(y1,y2):
    return np.sqrt(((y1 - y2)**2).mean())

def find_reg_val(f_train,f_val,y_train,y_val,soap,kernel,kmm_train,knm_train,sparse_train,regmin,regfac):
    x_val = soap.transform(f_val)
    reg     = regmin
    regs    = []
    scores  = []
    while True:
        regs.append(reg)
        # fit model
        weights = fit_gap_simple(f_train, kmm_train, y_train, knm_train, reg, jitter = 1e-40, target_type="Atom", solver="RKHS-QR")
        # predict on validation set and evaluate RMSE score
        y_pred = KRR(weights, kernel, sparse_train, units={'energy': 'ppm', 'length': 'AA'}).predict(x_val).flatten()
        scores.append(get_rmse(y_val,y_pred))
        # check if the minimum has been passed
        if len(scores) > 4 and np.all(scores[-1] > scores[-4:-1]):
            best_reg   = regs[np.argmin(scores)]
            best_RMSE  = np.min(scores)
            print("Best RMSE is :", best_RMSE," at a regularisation of :", best_reg)
            return best_reg
        # exponentially increase regularisation
        reg *= regfac

def train_val(f_train,f_val,y_train,y_val,kernel,kmm_train,knm_train,sparse_train,reg,n_models=8):
    x_val   = soap.transform(f_val)
    weights = []
    for j in range(n_models):
        # pick nrs CV training points from among the full training set and fit model
        ids = np.random.choice(np.arange(len(y_train)), int(len(y_train)/2), replace=False)
        weights.append(fit_gap_simple(f_train, kmm_train, y_train[ids], knm_train[ids], reg, jitter = 1e-40, target_type="Atom", solver="RKHS-QR"))
    # predict on validation set
    y_pred = KRR(weights,kernel,sparse_train,self_contributions=None,units={'energy': 'ppm', 'length': 'AA'}).predict(x_val).reshape((-1,len(weights))).T
    y_best = np.mean(y_pred,axis=0)
    y_err  = np.std(y_pred,axis=0)
    # evaluate rescaling factor for uncertainty estimates based on maximising likelihood of obtaining observed validation errors
    alpha  = -1.0/n_models + (n_models - 3)/(n_models - 1) * np.sqrt(np.mean((y_best - y_val)**2/y_err**2))
    # rescale uncertainty estimates
    y_err *= alpha
    return np.array(weights),alpha,y_best,y_err

# -----------------------------------------------
# model generation -- should not need to be touched
# -----------------------------------------------
# load and prepare training data
f_train = read('../Shielding_datasets/'+compound+'_ShiftML_frames_train_FPS_w_shifts.xyz', ':')
f_train,cs_train,nat_train = prep_frames(f_train,species)
print('loaded training frames')

# initialise features
hypers = {'soap_type' : "PowerSpectrum",
          'interaction_cutoff': rc,
          'max_radial': nmax,
          'max_angular': lmax,
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
hypers = get_optimal_radial_basis_hypers(hypers,f_train,expanded_max_radial=16)
print('initialised hypers')

# initalise SOAPs
soap,manager_train = calculate_features(f_train,hypers,auto_wrap=False)
print('initialised manager')

# load FPS order for environments
fps_indices = np.array([int(index) for index in np.loadtxt('../FPS_orders/env_fps_ids_'+compound+'_'+str(species)+'.dat')])
print('FPS-sorted')

# sparsify training points
sparse_train   = SparsePoints(soap)
sparse_indices = build_sparse_list(nat_train,fps_indices[:int(nat_train.sum()*fraction)])
sparse_train.extend(manager_train,sparse_indices)
print('sparsified')

# compute sparse and Gram kernel matrices
kern,kmm_train,knm_train = compute_kernels(soap,manager_train,sparse_train,soap_power=zeta,do_gradients=False,target_type="Atom")
print('kernels computed')

# load and prepare training data
f_val = read('../Shielding_datasets/'+compound+'_ShiftML_frames_val_w_shifts.xyz', ':')
f_val,cs_val,nat_val = prep_frames(f_val,species)
print('loaded validation frames')

# train models
opt_reg = find_reg_val(f_train,f_val,cs_train,cs_val,soap,kern,kmm_train,knm_train,sparse_train,regmin=1e-12,regfac=2)
weig,alpha,cs_best,cs_err = train_val(f_train,f_val,cs_train,cs_val,kern,kmm_train,knm_train,sparse_train,reg=opt_reg,n_models=nmodels)
print('models trained')

# save models
with open('../Shielding_models/'+compound+'_'+str(species)+'_rebuilt.pickle', 'wb') as f:
    pickle.dump([soap,kern,sparse_train,weig,alpha],f)
