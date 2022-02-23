#!/bin/bash

# initialize the socket and set up the simulation
~/i-pi/bin/i-pi input.xml > log.ipi &
sleep 30
~/i-pi/bin/i-pi-py_driver -u -a gapd-nvt-bI -m rascal -o ../../CSD_GAP_model.json,init.extxyz > log.gap &
~/dftb+/bin/dftb+ > log.dftb &
wait

