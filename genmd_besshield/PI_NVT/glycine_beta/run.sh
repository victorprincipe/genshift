#!/bin/bash

# initialize the socket and set up the simulation
~/i-pi/bin/i-pi input.xml > log.ipi &
sleep 10
~/i-pi/bin/i-pi-py_driver -u -a gapd-nvt-gb -m rascal -o ../../CSD_GAP_model.json,init.extxyz > log.gap &

~/dftb+/bin/dftb+ > log.dftb &
wait

