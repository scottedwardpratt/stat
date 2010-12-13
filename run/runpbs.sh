#!/bin/sh -login
#PBS -q main
#PBS -l mem=1000mb
set -e
echo irun=${PBS_ARRAYID};
./isHydro3 run${PBS_ARRAYID};
./b3d run${PBS_ARRAYID};


