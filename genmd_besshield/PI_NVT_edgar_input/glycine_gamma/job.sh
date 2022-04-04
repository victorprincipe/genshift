#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --partition=jobs
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks-per-node=64

module purge
module load intel
module load intel-mkl
module load intel-mpi

rm /tmp/ipi_dftb-nvt-gg /tmp/ipi_gapd-nvt-gg

i-pi input.xml > log.ipi &
sleep 10

mpirun -np 31 i-pi-py_driver -u -a gapd-nvt-gg -m rascal -o ../../CSD_GAP_model.json,init.extxyz > log.gap &

mpirun -np 31 ~/code/dftb+/bin/dftb+ > log.dftb &
wait
exit 0
