#!/bin/zsh
# ccs, cec24@phy.duke.edu
# run simple-test.R 
# currently hard-wired because i'm lazy
infile="../analysis/spectra-combined-s.dat"
errfile="../analysis/spectra-combined-errs-s.dat"
paramfile="../parameters/stats.param-combined-s"
nhold=5
export infile errfile paramfile nhold

R --slave --no-save < simple-test.R