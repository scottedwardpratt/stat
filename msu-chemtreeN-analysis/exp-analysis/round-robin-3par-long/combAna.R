## trains the emulator on the given round-robin folder.
## further functions can be used to generate the implausibility slices etc

## this is designed to be called by something that has defined mwString and
## is of course in the right place for the paths to work :)

source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings
initEmu()
source("~/local/include/emu-analysis/fnAnalysis.R")


## note that this script
## assumes that fileNames.R has already been loaded into the global scope
## fileNames.R defines: lumOutputFile, metOutputFile and so on

## load the model data
##
## first the luminosity data
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]
##
## now the metallicity data
modelDataMet <- abs(as.matrix(read.table(metOutputFile)))
nbinsMet <- dim(modelDataMet)[2]



if(nruns != dim(modelDataMet)[1]){
  stop("nruns modelDataMet doesn't match modelDataLum")
}

nbins <-  nbinsMet + nbinsLum 
modelData.big <- cbind(modelDataLum, modelDataMet)
modelData <- modelData.big


## load the design (design file is already in the global scope, from fileNames.R)
#designFile <- paste("./", mwString, "/design/design_", mwString.Up, ".dat", sep="")
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
## default point is given as last line, ignore it
designData.big <- as.matrix(read.table(designFile, col.names=desNames))#[1:nruns,]

designData <- designData.big


## load the experimental data, this is not what we actually want to do
## since we want to do a broad round-robin comparison we need to load the
## exp data for each of the runs we're going to generate comparisons against
## the exp-data is not used in the estimation process so we can start by setting this
## as a blank list
##
expData <- list()



rebuild <- 1
buffer <- paste("functional-data-", nruns, "-test-", mwString, ".dat", sep="")
if(rebuild == 1 || file.exists(buffer) == FALSE){
  ##
  ## generate a functional sample from the vars in global sope
  ## regression with a constant term seems to work better in general than
  ## a first order term
  fnData <- fn.sample.gen(cov.fn=1, reg.order=0)
  ## now do the pca decomp
  fnData <- fn.pca.gen(fnData, cutOff=0.98)
  ## estimate the thetas
  fnData <- fn.estimate(fnData)
  save(fnData, file=buffer)
} else {
  load(buffer)
}


