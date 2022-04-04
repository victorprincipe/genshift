#!/bin/bash
~/i-pi/bin/i-pi input.xml > log.ipi &
sleep 30
~/n2p2/src/interface/lammps-nnp/src/lmp_mpi < lmp_pbe_ts.in > log.lmp_pbe_ts &
wait
