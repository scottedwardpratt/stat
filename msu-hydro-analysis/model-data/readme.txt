madai-analysis/model-data

ccs, cec24@phy.duke.edu

The model-data folder is where the folders for each analysis are stored. Each analysis consists of a set of model runs and paremeters as well as data generated from this. 

=-=-=-=-=-=- PRE SETUP -=-=-=-=-=-=

An example analysis (./example) should contain the following folders, before running ./setup.sh

folders:
./example: 
./example/analysis - contains folders run1 ... run60 plus a default folder. Each run represents a location in the design space. Each run folder contains a set of centrailty folders, each of these in turn contains the results for that run and centrality in a results.dat file. 

./example/parameters - contains folders run1 ... run60 plus a default. Again each run reperesents a location in the design space. Each folder contains two files: stats.param and fixed.param. Stats.param contains the values that specify the location fo the run in the design (parameter) space and thus varies across runs. Fixed.param is the same in each run and lists the options which were fixed for the whole analysis.

=-=-=-=-=-=- POST SETUP -=-=-=-=-=-=


An example analysis with say 60 runs (called ./example) should contain the following after having run setup.sh. Setup will use the given specification templates and create combined output files, as well as creating a combined design.

The folders will stay the same but the following files should appear

files:

-- each type of file appears in cent classes 0to10, 10to20, 20to30, along with an errors file and a default file. Each file contains a header line listing the names of the variables and an additional column with the run-id.

spectra-combined-0to10.dat - the extracted star spectra from the results.dat files for all the runs in the first centrality class. 
spectra-combined-0to10-errors.dat - the errors for the same
spectra-combined-0to10-default.dat - the default "experimental" values to be used in comparisons

spec-and-flow-combined-0to10.dat - the extracted star spectra and v2 from results.dat

hbt-combined-0to10.dat - the extracted hbt data

./parameters/stats.param-combined-s - the combined and sorted design for the run. This should be used with any of the observables and any centrality.
./parameters/stats.param-default - the location of the default values. To be used when checking the exp-analysis



=-=-=-=-=-=- DURING ANALYSIS -=-=-=-=-=-=

During an emulator analysis the following files are created and stored here

files:

fn-data-<observable-name>-<cent-bin>-comb.dat - contains working data for the emulator, loading this in R will recover all the data needed to make predictions, this is an important file.



