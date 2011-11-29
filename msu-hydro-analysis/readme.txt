madai/madai-analysis 
C.Coleman-Smith, cec24@phy.duke.edu 29.11.2011

==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-==

Description:

Functions to train an emulator (estimate) and then generate predictions at arbitrary locations. Contains scripts and wrappers to support the MADAI MCMC code and to specifically integrate with the MSU 2d hydro code. Uses the emulator framework. 

==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-==

Depends: 

emulator project (supplies: libRbind emuRbind.R). Get the latest version from the github archive (https://github.com/jackdawjackdaw/emulator)
emu-analysis project (supplies: fnAnalysis.R). Get the latest version from the madai stats tree.

==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-==

Definitions of terms:

Each set of runs of the computer-model across a given design is referred to as an analysis. Each point in the design space is referred to as a run, the locations in the design space are often referred to as parameters. There are multiple centrality clasess per run. 

==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-==

Expected folder layout:

Each analysis is stored in a separate dir in ./model-data, each analysis dir should contain parameters and analysis directories. The analysis and parameters dirs each contain n run directories, run<1>...run<n>. Each run directory in analysis contains centrality class dirs, each of these contains a results.dat file giving the results for that centrailty.  Each run directory in parameters contains files which give the details of the location in the parameter space of that particular run.

==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-==

Contents:

data-prep - contains scripts and templates for extracting the analysis data, extraction of observables is controlled by 
template files (./data-prep/*-spec.dat) which can be easily modified. 

exp-analysis - folders and scripts for automatically generating plots of the observables against the design (plotplot.R) and  building an emulator and running implausibility analysis and making plots.

model-data - expected location for each analysis folder. See the readme and definitions above.

src - codes for training an emulator on a given analysis and for generating predictions.

setup.sh - extracts observables and design for a fresh analysis folder

==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-== ==-=-=-==

Getting Started:

This is just enough information to get an MSU analysis running.

0) extract the analysis.tar and parameters.tar folders in a common folder. eg ./model-data/analysis-test

1) Processing data for analysis:

Run setup.sh with the path to the folder you created in step 0. This should extract analysis and parameters files

eg: 
./setup.sh ./model-data/analysis-test/

See data-prep/readme.txt for more detailed instructions on extracting the data.

2) Building the emulator model:

src/makeThetas.sh will generate hyper-params for the data. Arguments are: 1 path to the analysis folder created in step 1, path to a model data file extracted by setup.sh in step 1. I.e hbt-combined-0to10.dat etc.

eg:

./src/makeThetas.sh ./model-data/analysis-fall-test/ ./model-data/analysis-fall-test/hbt-combined-0to10.dat 


3) Predicting locations

src/computePoints.sh will make predictions at a list of points. 

./src/makeThetas.sh ./model-data/analysis-fall-test ./model-data/analysis-fall-test/hbt-combined-0to10.dat




