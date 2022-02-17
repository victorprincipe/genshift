#!/usr/bin/env python

# -----------------------------------------------
# imports
# -----------------------------------------------
import numpy as np
from ase.io import read, write
import pickle
import rascal
from rascal.neighbourlist.structure_manager import mask_center_atoms_by_species
from rascal.models import KRR

# -----------------------------------------------
# specify inputs variables
# -----------------------------------------------
species  = [1,6] # list of chemical species for which to compute shieldings
compound = 'benzene' # compound, required for loading the correct ML shielding model
nmodels  = 16 # all committees are composed of 16 individual models
inname   = 'no_eq_simulation_ase.xyz' # filename/path of the input trajectory to calculate shieldings for
outname  = 'no_eq_simulation_ase_w_cs.xyz' # filename/path of the file to which the trajactory with shieldings should be written
nbatch   = 500 # number of environments to be considered at a given time (to limit memory requirements)

# -----------------------------------------------
# function definitions -- should not need to be touched
# -----------------------------------------------
def predict(f_test,soap,kern,feat,weights,alpha):
    x_test   = soap.transform(f_test)
    print('features evaluated')
    y_pred   = KRR(weights,kern,feat,self_contributions=None,units={'energy': 'ppm', 'length': 'AA'}).predict(x_test).reshape((-1,len(weights))).T
    print('models evaluated')
    # rescaling of differences and corresponding uncertainty estimates according to likelihood maximisation for obtained observed validation errors
    y_final  = y_pred * 0.0 + np.mean(y_pred,axis=0)
    y_final += alpha * (y_pred - np.mean(y_pred,axis=0))
    return y_final

# -----------------------------------------------
# shielding calculation -- should not need to be touched
# -----------------------------------------------
# load up test frames
traj = read(inname,':')
# prepare test frames
bohr2ang = 0.529177
for ifrm,frm in enumerate(traj):
    # set arrays in which the chemical shieldings, uncertainty estimates, and individual predictions from the members of the committee will be stored
    frm.set_array('CS',np.zeros(frm.get_global_number_of_atoms()))
    frm.set_array('CSerr',np.zeros(frm.get_global_number_of_atoms()))
    frm.set_array('CSens',np.zeros((frm.get_global_number_of_atoms(),nmodels)))
    # convert cell parameters and atomic positions to Angstrom
    frm.positions *= bohr2ang
    frm.cell      *= bohr2ang
    # wrap all atoms into the unit cell
    frm.wrap(eps=1e-11)
print('frames loaded and prepared')
# run over chemical species
for sp in species:
    print('considering chemical species',sp)
    # load ShiftML model
    with open('/local/home/local/Documents/principev/edgar_reconstruction/Shielding_models/'+compound+'_'+str(sp)+'.pickle','rb') as f:
        [soap,kern,feat,weig,alpha] = pickle.load(f)
    print('loaded shielding model')
    # mask all atomic centers except for species at hand
    nsp = []
    for ifrm,frm in enumerate(traj):
        frm.set_array('center_atoms_mask',np.zeros(frm.get_global_number_of_atoms(), dtype=bool))
        mask_center_atoms_by_species(frm,species_select=[sp,])
        nsp.append(np.sum(frm.get_array('center_atoms_mask')))
    # split trajectory into batches
    ibatch = [0]
    jbatch = []
    for ib in range(int(np.sum(nsp)/nbatch + 0.5)):
        ifrm = 0
        msp = 0
        while msp < (ib+1) * nbatch:
            msp = np.sum(nsp[:ifrm])
            ifrm += 1
            if ifrm == len(traj):
                ifrm += 1
                break
        jbatch.append(ifrm)
        ibatch.append(ifrm)
    print('split frames into',len(jbatch),'batches')
    # predict for test frames
    for ib in range(int(np.sum(nsp)/nbatch + 0.5)):
        pred      = predict(traj[ibatch[ib]:jbatch[ib]],soap,kern,feat,weig,alpha)
        print('predicted for batch',ib+1)
        pred_mean = pred.mean(axis=0)
        pred_err  = pred.std(axis=0)
        # save trajectory with predicted CS for future use
        counter = {}
        counter[sp] = 0
        for ifrm,frm in enumerate(traj[ibatch[ib]:jbatch[ib]]):
            # get CS
            cs    = frm.get_array('CS')
            cserr = frm.get_array('CSerr')
            csens = frm.get_array('CSens')
            # set predicted CS for H and C
            mask         = np.where(frm.get_atomic_numbers() == sp)[0]
            cs[mask]     = pred_mean[counter[sp]:counter[sp]+len(mask)]
            cserr[mask]  = pred_err[counter[sp]:counter[sp]+len(mask)]
            csens[mask]  = pred.T[counter[sp]:counter[sp]+len(mask)]
            counter[sp] += len(mask)
            # set collected CS for frame
            frm.set_array('CS',cs)
            frm.set_array('CSerr',cserr)
            frm.set_array('CSens',csens)
        print('predictions stored for batch',ib+1)
    write(outname,traj)
