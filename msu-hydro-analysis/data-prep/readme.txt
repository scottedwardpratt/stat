data-prep/readme.txt
ccs cec24@phy.duke.edu 2011-March/April (MSU)

This folder contains scripts needed to prepare output from the current RHIC analysis for use with the emulator. They are used by setup.sh, you shouldn't need to ever call them directly.

To extract different observables you should make a copy of everything-spec.dat and remove (or comment with #'s) out lines you don't want. Then follow the steps in setup.sh to add a new extraction, or run generate-analysis-files.rb directly.

files:
generate-analysis-files.rb - extracts a given set of observables (supplied on stdin) from an analysis folder
process-parameters.sh - extracts the combined design from a given parameters folder
testSpec.dat - example observable extraction template 
spectra-only-spec.dat - observable extraction file, picks out all the star spectra observables
star-spectra-and-flow.dat - observable extraction file, picks out all the star spectra and the v2 measurements
hbt-only-spec.dat - observable extraction file, picks out all the star hbt data.
everything-spec.dat - observable extraction file, picks out all possible fields, useful as a starting point, delete out rows you don't want for your own custom specificationp


folders:
proc-param - contains the helper scripts run by process-parameters to extract the combined design


