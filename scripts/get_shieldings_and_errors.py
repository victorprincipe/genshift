import sys
import ase.io
import numpy as np
import pandas as pd

xyz = sys.argv[1]
save = sys.argv[2]

def get_shieldings_and_errors(xyz_file, save_filename):
    '''
    Usage: python get_shieldings_and_errors.py (xyz_file) (save_filename)
    Read an ASE-formatted xyz trajectory, get the chemical shieldings and errors, and save them to a new CSV file.
    The file will be saved as save_filename.csv (i.e. do NOT add .csv to the save_filename option)
    '''
    
    #get number of frames in traj
    all_frames= len(ase.io.read(xyz_file, index=':'))
    print('Frames to read: '+ str(all_frames))
    
    #get array of chemical symbols for columns
    symbols=ase.io.read(xyz_file, index=1).get_chemical_symbols()
    
    #number atoms of each species
    c = 0
    h = 0
    n = 0
    o = 0
    for num in range(len(symbols)):
        if symbols[num] == 'C':
            symbols[num] = symbols[num]+str(c+1)
            c += 1
        elif symbols[num] == 'H':
            symbols[num] = symbols[num]+str(h+1)
            h += 1
        elif symbols[num] == 'N':
            symbols[num] = symbols[num]+str(n+1)
            n += 1
        elif symbols[num] == 'O':
            symbols[num] = symbols[num]+str(o+1)
            o += 1
    
    #create array of column headers for errors
    error_symbols = []
    for symbol in symbols:
        error_symbols.append(symbol + 'err')
    
    print('Atoms per frame: '+ str(len(symbols)))
    
    shieldings_raw=[]
    shielding_errors_raw=[]
    
    atoms=ase.io.read(xyz_file, index=':')
    
    #get chemical shieldings and errors form each frame
    if atoms[0].has('cs_iso')==True:   #check if array is indexed using Matthias's method
        for i in range(len(atoms)):
            cs=atoms[i].get_array('cs_iso')
            shieldings_raw.append(cs)
    else:
        for i in range(len(atoms)):
            cs = atoms[i].get_array('CS')
            cserr = atoms[i].get_array('CSerr')
            shieldings_raw.append(cs)
            shielding_errors_raw.append(cserr)

    #create numpy arrays of shieldings and errors
    shieldings=np.array(shieldings_raw)
    shielding_errors=np.array(shielding_errors_raw)
    
    #save shieldings and errors to new CSV
    df_shieldings = pd.DataFrame(shieldings, columns = symbols)
    df_total = df_shieldings
    
    if len(shielding_errors)>0:
        df_errors = pd.DataFrame(shielding_errors, columns = error_symbols)
        for column in df_errors.columns:
            df_total[column] = df_errors[column]
    
    df_total.to_csv(save_filename+'.csv')
    print('Shieldings and errors saved as: '+save_filename+'.csv')
    
    return

get_shieldings_and_errors(xyz, save)
