import sys
import ase.io

'''
Usage: python discard_equilibration.py (.xyz file to read, percentage = 5) 
Functionality: Discard the equilibration period of a simulation. This is equal to the first (percentage) number of frames. If no value is passed for the percentage, the first 5% will be discarded by default. The trajectory must be in extended ASE format. The new file is saved in the same format with 'no_eq_' in front of the input file name.
'''

#set default percentage of frames to discard if no second argument is passed
if len(sys.argv)==2:
    xyz= sys.argv[1]
    percent = 5

# set the passed percentage of frames to discard
elif len(sys.argv)==3:
    xyz= sys.argv[1]
    percent= int(sys.argv[2])
else:
    print("Incorrect number of arguments passed (must be 1 or 2). Check 'Usage' section.")

def discard_equilibration(xyz_file, percentage):

    #Read xyz file
    sim = ase.io.read(xyz_file, index = ':')
    print('Total number of frames: '+ str(len(sim)))
    
    #Index frames to chuck
    to_chuck= int(len(sim)*(percentage/100))
    print('Discarding first ' + str(percentage) + '% of frames.')
    
    #Remove unwanted frames
    to_use= sim[to_chuck:]
    print('Frames remaining:'+ str(len(to_use)))

    #Write xyz file (second argument passed)
    ase.io.write('no_eq_'+str(xyz), to_use)

    print('File saved as:', 'no_eq_'+str(xyz))
          
    return

discard_equilibration(xyz, percent)
