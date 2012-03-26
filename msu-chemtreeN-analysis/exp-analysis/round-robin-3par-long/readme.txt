21.02.2012
ccs, cec24@phy.duke.edu

steps to provision a new folder

1) copy in ../round-robin/*.R 
2) copy in ../round-robin/*.rb
3) mkdir ./raw-data
4) copy in raw data
5) edit setup.rb
	 - need to insert correct runName suffix, in this case '_3par.dat'
	 - make sure that the designName suffix is right too, in this case '_3par_sorted.dat'
6) run setup.rb -> this creates the mw-round folders and populates their subfolders with design and obesrvables
7) run setup-r-stubs.rb -> this creates r-stubs to build the emulator, make implausibility grids etc
8) edit fileNames.R 
	 - ensure that the _3par _5par or whatever filename suffix is correct

building the emulator

1) run ./rpt-combAna.sh -> this will train emulators for each mwI-round data set
2) run ./rpt-implaus.sh -> this will create implausibily data and output implausibility grids, this will				
	 make nice implausibility plots
3) run ./post-plot.sh -> this creates posterior plots, Exp[-I^2 / 2], and dump an evidence matrix
	 
