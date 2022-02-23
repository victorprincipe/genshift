#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --partition=jobs
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=64


rm /tmp/ipi_dftb_nvt /tmp/ipi_gap_d_nvt

i-pi input.xml > log.ipi &
sleep 30

mpirun -np 32 i-pi-py_driver -u -a gapd-nvt-bI -m rascal -o ../../CSD_GAP_model.json,init.extxyz > log.gap &

mpirun -np 32 ~/code/dftb+/bin/dftb+ > log.dftb &
wait
