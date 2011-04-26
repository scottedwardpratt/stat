rhic-data-prepare/readme.txt
ccs cec24@phy.duke.edu 2011-March/April (MSU)

This folder contains scripts needed to prepare output from the current RHIC analysis for use with the emulator.
Data for emulator use is assumed to be in a subfolder of ${MADAI_ROOT}/model-data, this file should be in 
${MADAI_ROOT}/src/rhic-data-prep, on my laptop ${MADAI_ROOT} is ~/projects/madai-analysis/.

The rhic-data for a given run should consist of analysis and parameters tarball. Unzip these into a subfolder of 
the ${MADAI_ROOT}/model-data folder. 

----the analysis folder-----
Contains N_runs subfolders labelled run_1...run_N. Each run_i subfolder contains N_qualifiers subfolders again, in the RHIC case the qualifiers are the impact parameters for each run. There are 5 impact parameters. Within the qualifier subfolders is a results.dat file, this contains a list of the output values of the model in the format
<value_type> <value_name> <value_value> <value_error>.

To get this data processed into a usable format,Run the wrapper process-analysis.sh with the first argument as the path to the top of the analysis file-tree that you want to process. This will create sorted files called: 
spectra-combined-s.dat, 
spectra-combined-errs-s.dat,
flow-combined-s.dat,
flow-combined-errs-s.dat,
hbt-combined-s.dat,
hbt-combined-errs-s.dat.

Each line in the file contains a subset of the results for ALL qualifiers for a given run. The spectra-combined files only reports the first 3 moments of the spectra, the flow-combined file reports the first 2 moments, etc. The errors files follow the same format, but contain the approximate prior error on the corresponding values.

These can now be fed into ${MADAI_ROOT}/src/makeThetas.sh to begin the estimation.

--- the parameters folder----
Contains N_subfolders, again run_1...run_N. For the emulator process we only care about the stats.param file within each folder. 

These fiels are all collated together, one line for each run, into the file stats.param-combined-s in the top of the paramters folder.

Run process-parameters.sh with $1 as the path to the parameters folder to produce this.



