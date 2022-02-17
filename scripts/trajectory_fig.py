import sys
import ase.io
import pandas as pd
import matplotlib.pyplot as plt

'''
Usage: python trajectory_fig.py (simulation.out file) (percentage of frames to use)
Takes a simulation.out file and produces plots of various properties over the course of the trajectory, saved as 'traj_properties.png'.
You can pass a second argument, which is equal to the percentage of frames (from first frame) that will be used for the plot (default = 100)
NOTE: it is assumed that the simulation.out file contains the columns in the 'header' variable in the same order.
If this isn't the case, please change the headers line above, and also make sure that the usecols option
below reflects the number of columns that are present.
'''

analyse = sys.argv[1]

if len(sys.argv)==3:
    percentage = int(sys.argv[2])
else:
    percentage=100

def trajectory_fig(sim_outfile, percent):
    
    # Set column headers
    headers=['Step', 'Time (ps)', 'Pot_component_raw(0)','Pot_component_raw(1)', 'Temperature', 'Conserved']
    
    # Read simulation file
    sim= pd.read_csv(sim_outfile, sep=' ', skipinitialspace=True, comment='#', usecols=[0,1,2,3,4,5], names=headers)
    
    # Index frames to use for plot 
    frames=int(percent*len(sim)/100)

    #Create plots
    plotted=sim.iloc[0:frames].plot(kind='line', subplots=True, x='Step', figsize=(10,20))
    fig=plotted[0].get_figure()
    
    #Save plot
    fig.savefig('traj_properties.png')
    
    print('Figure saved as traj_properties.png')
    
    return
 
trajectory_fig(analyse, percentage)
