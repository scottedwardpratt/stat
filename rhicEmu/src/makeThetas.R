## makeThetas:
## ccs, cec24@phy.duke.edu 2011-March/April (msu)
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
## errfile -> same as infile but contiaining the estimated errors for each training point (not used yet)##
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
source("gen-samp.R") # load the sample set
source("emu-pca.R") # load the pca fns
# note, need to change this to .so or whatever for deployment
dyn.load("~/local/lib/libRBIND.dylib") #
# move these files into the project root is  probably best
source("~/local/include/libRbind/EmuRbind.R")
source("~/local/include/libRbind/multivar.R")


# this file reads the given sample set and generates a pca.decomp for it and the thetas

cargs <- Sys.getenv(c('infile', 'errfile', 'paramfile', 'outpath'))
inSamples <- cargs[1]
inSamplesErrs <- cargs[2]
inParams <- cargs[3]
outPath <- cargs[4] 

## might want to print some blurb here?
## like this is bloopy-bloop, herping and derping
cat("#Estimating Thetas \n")
cat("#OutPath ", outPath, "\n")
cat("#Intermediate data will be written to: SampleTemp.dat, PcaTemp.dat\n")
cat("#Final Hyper-Params are written to: theta-table.dat\n")

des <- load.design(inParams)
sample <- load.sample(inSamples, inSamplesErrs, center=FALSE, centerWithErrs=TRUE)

# save the full sample for later use
buffer <- paste(outPath, "/SampleTemp.dat", sep="")
save(sample, file=buffer)

# now construct the sample set that the pca wants
tvec <- seq(1, sample$ntps)
samplePCA <- list(t=tvec, y=t(sample$y), des=t(des$des))

# do the pca decomp
pca.decomp <- emu.pca(samplePCA, nr=4, silent=TRUE)
buffer <- paste(outPath, "/PcaTemp.dat", sep="")
save(pca.decomp, file=buffer)

# finally estimate the thetas
runTime <- system.time(thetas.est <- estimateThetas(pca.decomp))
buffer <- paste(outPath, "/theta-table.dat", sep="")
write.table(thetas.est, buffer, row.names=FALSE, col.names=FALSE) # and dump to file

cat("#RunTime ", runTime, "\n")

