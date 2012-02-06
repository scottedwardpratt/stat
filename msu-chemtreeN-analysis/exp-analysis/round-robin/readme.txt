each folder at the root contains a functional-data-<mwname>.dat file created by <mwName>CombAna.R

this is then used to make comparisons against the other mw observations in each of the mw1, mw2, etc sub folders

to recreate everything
1) run setup-folder.rb - copies data from raw-data, exp-output into the right mw folders
	 filenames have to be exactly the same
2) run rpt-combAna.sh - trains the emulator, this works by running the script combAna.R which sets up the model and		calls all the R-fns
3) run rpt-implaus.sh - creates implausibility plots
