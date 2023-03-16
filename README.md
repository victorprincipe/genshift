# Genshift - Predicting Chemical Shifts from ML-Enhanced MD simulations
Here you will find:
-
-
-
-

## Chemical Shieldings prediction workflow

### Producing input files from downloaded COD structure (.cif) files
- The Crystallography Open Database (COD) has lots of structures with downloadable .cif files
- The .cif files need to be converted to .xyz and .cfg files, which can be done using the ASE package
- These .cif files can in turn be converted to LAMMPS input files (.lmp ) using the command atomsk (filename).cfg lammps
	- The LAMMPS input files are saved as init.data in Edgar's PIMD_examples folders! But they are actually the same format as .lmp
- When starting simulations, LAMMPS gives the following warning:
```bash
CAUTION: Please ensure that this mapping between LAMMPS
atom types and NNP elements is consistent:
---------------------------
LAMMPS type | NNP element
---------------------------
1 <-> H ( 1)
2 <-> C ( 6)
3 <-> N ( 7)
4 <-> O ( 8)
---------------------------
```
- Make sure that the numbering of atom types in the .lmp file is the same as what is shown above- they are usually different when the .lmp file is produced!!
- You also need to change the init.xyz files
	- Manually add the following line with lengths and angles
```
# CELL(abcABC): 5.223000 12.435000 5.18862100 90.00000 111.14000 90.00000 Step: 45 Bead: 0 positions{angstrom} cell{angstrom}
```
- (And delete this part from the init.xyz file as it's not necessary):
```
Lattice="5.223 0.0 0.0 0.0 12.435 0.0 -2.00628507308967 0.0 5.1886211275730645" Properties=species:S:1:pos:R:3 spacegroup="P 21/c" unit_cell=conventional occupancy="_JSON {\"0\": {\"N\": 1}, \"1\": {\"H\": 1}, \"2\": {\"H\": 1}, \"3\": {\"H\": 1}, \"4\": {\"C\": 1}, \"5\": {\"H\": 1}, \"6\": {\"H\": 1}, \"7\": {\"C\": 1}, \"8\": {\"O\": 1}, \"9\": {\"O\": 1}}" pbc="T T T"
```

### Workflow for Bespoke MD + Bespoke Shielding:

1. Run PIMD NVT/NpT simulations
2. Use the `trajectory_fig.py script` on the simulation.out file in order to determine the equilibration period, and to see if the simulation was stable (look at the conserved energy quantity)
3. Convert a pos_*.xyz  file to ASE format using the provided `convert_ipi_to_ase.py`  script
4. Use the `discard_equilibration.py`  script to discard the equilibration frames from the ASE-converted pos_*.xyz  file
5. Calculate chemical shieldings using the `calculate_shieldings.py`  script
6. Use the `get_shieldings_and_errors.py` script to retrieve data of shieldings and errors from the xyz file containing shielding predictions
7. Use the `plottings.ipynb` notebook to process results and plot graphs
	- You might need to remove rows with zeros from the dataframe (read in through the .csv file) using the remove_zeros function


### Running Bespoke MD simulations using Lammps + n2p2:

- Don't compile Lammps- as n2p2 has an option for compiling Lammps separately using  make lammps-nnp in the src directory
- However, Lammps requires the i-PI socket interface to work with i-PI, which is part of the optional user-misc package in lammps
- So, after doing make lammps-nnp , LAMMPS is installed in the interface directory of n2p2 and unpacked into lammps-nnp
- In the lammps-nnp src directory, you can make yes-user-misc to activate the installation of the user-misc package
- Then, you can do a simple make mpi  to compile lammps, where the binary lmp_mpi will be added to the either the src or bin directory of lammps-nnp (don't remember which one it's copied to)
- Then, lastly, you have to delete the copy of lmp_mpi from the n2p2/bin directory, and then copy the new lmp_mpi from the src directory of lammps-nnp to n2p2/bin


### Running General MD simulations using Aditi's potential:

- Requires installation of DFTB+ with the following options set to TRUE in the cmake.configure file EXCEPT WITH_MPI should be set to FALSE:

![dftb_compile](https://user-images.githubusercontent.com/75172693/225672601-17f998dd-3666-428f-901c-cc13ecbece9f.png)

- However, if you're compiling for use on the cluster, you should use WITH_OMP set to FALSE as well
	- When compiling on the cluster, make sure to download the optional externals using ./utils/get_opt_externals
- For the DFT in file (dftb_in.hsd):
	- You have to add/remove hubbard derivs and angular momenta that aren't relevant to the molecule (e.g. no N or O in benzene)
	- The input file requires specific k-points for each compound, which can be calculated using the script in compute_kpoints.ipynb in genshift/genmd_besshield 
	- Only the first three numbers for k_points need to be changed (not the weights, e.g. 1.0)
- The librascal potential needs the init.extxyz format -- look at look at compute_kpoints.ipynb for the function to compute it 
- init.gen  can be created from init.xyz using the ASE package 
	- You might need to read the init.xyz using the function read_xyz as written by Emma in compute_kpoints.ipynb !!
- You might have to pip install mkl  before running sims
- Make sure to use the 500 MB CSD_GAP_model.json ! (It is only 23KB on Emma's github)
- You should change the input.xml files to:
	- Reduce the number of nbeads - this will have a big effect on simulation efficiency, but may reduce the accuracy (Emma used a value of 1...)
		- Needs to be >1, otherwise you're just doing classical sampling, not PIMD!
	- Make sure that the two ffsocket name have "driver-dftb" and "driver'gap"
	- Change forces to have the two force forcefield options equal "driver-dftb" and "driver-gap"
	- Optional, but change the properties to have step, time{picosecond}, potential{electronvolt}, conserved, temperature{kelvin}, ensemble_temperature{kelvin}
- Make sure that socket driver addresses match in input.xml , dftb_in.hsd ,  and run.sh



### Using Matthias's general shielding potential on the bespoke simulations:

- The ASE-converted xyz trajectory files, with no equilibration period, were copied from the bespoke MD + bespoke shielding directories
	- Look in these original folders for any raw files that might be necessary, as wel as the figure of the trajectory properties over time
- The ASE-converted xyz files from the bespoke MD simulations first have to be converted (from Bohr) to units of angstroms!
	- This can be done using the `bohr_2_ang.py` script 
- Use the script `calculate_shieldings_general.py`  to use Matthias's model for shielding prediction
	- Once you apply Matthias's model, the shielding predictions are stored in an array called cs_iso
	- Results can be analysed using the `plottings.ipynb` in the home directory of genshift
- There isn't any prediction of the error from the ML method like Edgar's model has.



### Plotting results:

- Look at plottings.ipynb in the `plots/` directory of genshift
	1. Read chemical shifts from the .csv file with chemical shieldings for each nucleus 
	2. If needed, remove empty rows from this dataframe using the function remove_zeros
	3. Create lists of all chemical shielding values from chemically equivalent nuclei and add them to a new dataframe
		- Can plot histograms of these results with  `test_processed.hist(column='C1', bins=200)`
		- Can also plot a density line graph (as done in the paper) by estimating the density function using the gaussian_kde() method from scipy.stats , as shown below:

```python
normalized_c1 = gaussian_kde(test_processed['C1'])
xs = np.linspace(min(test_processed['C1']),max(test_processed['C1']),200)
normalized_c1.covariance_factor = lambda : .25
normalized_c1._compute_covariance()
plt.plot(xs, normalized_c1(xs))
```

- Where 'C1' is the list (or array) of benzene C1 chemical shieldings as part of the test_processed dataframe
	
	4. Find mean, stddev and errors
	5. Convert from chemical shielding to chemical shift
		- This can only be done if there are multiple chemical environments for a species. In the case of N, as there is only one environment, shielding=shift.


Use the following commands once simulations have been run:

To clean up files and folders (only necessary if there are a lot of replicas):

```bash
mkdir sim_frc
mkdir sim_pos
mkdir sim_restart
mv simulation.frc_* sim_frc/
mv simulation.pos_* sim_pos/
mv simulation.restart_* sim_restart/
cp sim_pos/simulation.pos_00.xyz .
```

Then, to obtain a figure of the trajectory properties, convert simulations to ASE format, and discard equilibrations:

```
python ~/projects/genshift/scripts/trajectory_fig.py simulation.out
python ~/projects/genshift/scripts/convert_ipi_to_ase.py simulation.pos_00.xyz simulation_ase.xyz
python ~/projects/genshift/scripts/discard_equilibration.py simulation_ase.xyz
```

1. If using the bespoke shielding models, issue the following commands:
```
nohup python ~/projects/genshift/scripts/calculate_shieldings_XXXX.py &
```
- Where XXXX is either "benzene", "glycine", or "succinic"

2. If using Matthias's model, issue the following commands:
```
python ~/projects/genshift/scripts/bohr_2_ang.py no_eq_simulation_ase.xyz
python ~/projects/genshift/scripts/calculate_shieldings_general.py no_eq_simulation_ase_ang.xyz
```

Then, to obtain only the shieldings in a csv, use:
```
python ~/projects/genshift/scripts/get_shieldings_and_errors.py no_eq_simulation_ase_w_cs.xyz shieldings
```
