#!/bin/sh -login
#PBS -q main
#PBS -l mem=1000mb
#PBS -l walltime=1:00:00
set -e
echo irun=${PBS_ARRAYID};
cd /mnt/home/prattsc/git/rhic/run;
./isHydro3 run${PBS_ARRAYID};
./b3d run${PBS_ARRAYID};


