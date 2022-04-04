import sys
import ase.io
import pandas as pd
import matplotlib.pyplot as plt

'''
Usage: python trajectory_fig.py (simulation.out file) (percentage of frames to use)
Takes a simulation.out file and produces plots of various properties over the course of the trajectory, saved as 'traj_properties.png'.
You can pass a second argument, which is equal to the percentage of frames (from first frame) that will be used for the plot (default = 100)
NOTE: please make sure that the usecols option below, during reading of simulaiton file, reflects the number of columns that are present.
'''

analyse = sys.argv[1]

if len(sys.argv)==3:
    percentage = int(sys.argv[2])
else:
    percentage=100

def trajectory_fig(sim_outfile, percent):
    
    # Get column headers
    headers = []
    with  open(sim_outfile,'r') as cmt_file:    # open file
        for line in cmt_file:    # read each line
            if line[0] == '#':    # check the first character
                line = line[21:]    # remove first '#'
                para = line.split(' :')     # seperate string by ':'
                if len(para) == 2:
                    headers.append(para[0])
    
    # Read simulation file
    sim= pd.read_csv(sim_outfile, sep=' ', skipinitialspace=True, comment='#', usecols=[0,1,2,3,4,5], names=headers)
    
    # Index frames to use for plot 
    frames=int(percent*len(sim)/100)

    #Create plots
    plotted=sim.iloc[0:frames].plot(kind='line', subplots=True, x='step', figsize=(10,20))
    fig=plotted[0].get_figure()
    
    #Save plot
    fig.savefig('traj_properties.png')
    
    print('Figure saved as traj_properties.png')
    
    return
 
trajectory_fig(analyse, percentage)
