## makeThetas:
## ccs, cec24@phy.duke.edu 2011-March/April (msu)
## heavy modifications made in 2011-October (msu)
## switched from using emu-pca.R to using fnAnalysis.R which is now installed by the
## emu-analysis project
##
## reads analysis and parameter data, does a pca decomposition and produces a set of  estimated
## hyperparameters (thetas) which can be used to emulate the given data (by computePoints)
##
## env: infile, errfile and paramfile, strings which give the locations of the
## input, error and parameter files to used
##
## infile -> an ASCII file (ordered by run id) containing training values (simulation generated data) 
## for each set of parameters. Currently there is a header which explains the columns for 1 impact 
## For other applications this could be simpler / harder
##
## paramfile -> an ASCII file (ordered by run id, or at least in the same was as infile) containing
## the locations in paramspace used to generate each of the corresponding training poitns
##
## outpath -> the location where we will save the thetas, pca decomp and samples for further use
##
## (lines in infile)
##
## output: saved pca.decomp R-structure (see emu-pca.R for more details)
##         saved sample R-structure (see gen.sample for more details)
##         ascii-table of optimum hyperparameters theta-table.dat
##
## the hyperparameters are saved as ascii as it's sometimes interesting to actually look at them or
## one might want to actually use them in some custom analysis. The R-data files are really only
## useful for computePoints
##
## 
#

## not needed?
##source("gen-samp.R") # load the sample set

source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings
# now load the emulator lib
#initEmu()

arch <- system("uname -s", intern=TRUE)
if(is.loaded("callEstimate") == FALSE){
  libNAME <- "~/local/lib/libRBIND"
  if(arch == "Linux"){
    libNAME <- paste(libNAME, ".so", sep="")
  } else if(arch =="Darwin"){
    libNAME <- paste(libNAME, ".dylib", sep="")
  } else {
    buffer <- paste("error: uname -s gives ", arch , " not supported", sep="")
    stop(buffer)
  }
  dyn.load(libNAME)
}


## load the functional analysis routines
source("~/local/include/emu-analysis/fnAnalysis.R")

## read the input params from the shell script
cargs <- Sys.getenv(c('infile', 'paramfile', 'outpath'))
modelOutput.file <- cargs[1]
designFile <- cargs[2]
outPath <- cargs[3]


analysis.name <- basename(modelOutput.file) ## this should be set automatically from infile
save.buffer <- paste(outPath, "/fn-data-", analysis.name, sep="")

cat("#Estimating Thetas \n")
cat("#OutPath ", outPath, "\n")

cat("#full data will be written to: ", save.buffer, "\n")

## load the model (strip off the first col which is the run id)
modelData <- as.matrix(read.table(modelOutput.file, header=TRUE))

## the number of reps of the model we have 
nruns <- dim(modelData)[1]
## the number of observables 
nbins <- dim(modelData)[2]
cat("#nruns: ", nruns , " nbins: ", nbins, "\n")

## load the design, need to re-label as they change
desNames <- c("GLAUBER_WNBIN_RATIO", "GLAUBER_K_TAU",
              "HYDRO_INIT_NS", "HYDRO_INIT_FLOW", "HYDRO_SVRATIO",
              "EQOFST_BUMP_HEIGHT")

## again strip off the run id
designData <- as.matrix(read.table(designFile)[,-1])
dimnames(designData) <- list(rownames=NULL, colnames=desNames)

## we would normally load the experimental data at this stage, but since
## the purpose of this program is only to generate thetas we don't need to
## the default run serves as our comparison to start with...

expData <- c()

## generate the fn sample from our global vars
fnData <- fn.sample.gen(cov.fn=1, reg.order=0)
## do the pca.decomp
fnData <- fn.pca.gen(fnData, cutOff=0.95)
## estimate the thetas
runTime <- system.time(fnData <- fn.estimate(fnData))
cat("#RunTime ", runTime, "\n")  
save(fnData, file=save.buffer)
buffer <- paste(strsplit(analysis.name, ".", fixed=TRUE)[[1]][1], "-thetas.txt", sep="")
write.table(fnData$thetas.est, buffer, row.names=FALSE, col.names=FALSE) # and dump to file




