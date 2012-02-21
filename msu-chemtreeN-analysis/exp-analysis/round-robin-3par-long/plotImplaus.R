## use the trained emulator (mwActual0 data to generate implausibility plots
## for the observations of mwCompare

source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings
initEmu()
source("~/local/include/emu-analysis/fnAnalysis.R")

## 21.02.2012
## again this code assumes that fileNames.R has been loaded,
## fileNames defines: lumOutputFile, metOutputFile designFile, expDataFileLum, and expDataFileMet


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

nbins <-  nbinsMet + nbinsLum 
modelData <- cbind(modelDataLum, modelDataMet)

## load the design
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
## default point is given as last line, ignore it
designData <- as.matrix(read.table(designFile, col.names=desNames))[1:nruns,]

## now load the observations of the galaxy we want to compare AGAINST
#expDataFileLum <- paste("./", mwStringCompare, "/lum_fun_observed_", mwStringCompare.Up, ".dat", sep="")
#expDataFileMet <- paste("./", mwStringCompare, "/metallicity_MV_observed_", mwStringCompare.Up, ".dat", sep="")

expDataLum <- as.matrix(read.table(expDataFileLum))
expDataMet <- as.matrix(read.table(expDataFileMet))

expData <- list(obsValue=c(expDataLum[2,], abs(expDataMet[1,])),
                obsError=c(expDataLum[3,], abs(expDataMet[2,])))

## now we load the fn data
buffer <- paste("functional-data-", nruns, "-test-", mwString, ".dat", sep="")
load(buffer)

## now wodge in the expData, but actually the fn fn.implaus.slice uses the globally available expData
fnData$exp.obs = expData

## systematically generate all projections of the joint implaus,
## this will make sections through the implaus for any number of parameters (neat)
combs <- combn(nparams, 3)
ncombs <- dim(combs)[2]
for(i in 1:ncombs){
  fixedValVec <- rep(0, nparams)
  fixedValVec[combs[1,i]] <- NA
  fixedValVec[combs[2,i]] <- NA
  fixedValVec[combs[3,i]] <- NA
  stepData <- fn.emu.steps(fnData, combs[1,i], combs[2,i], combs[3,i], fixedValVec)
  impSteps <- fn.implaus.steps(fnData, stepData)
  buffer <- paste("./", mwStringCompare, "/images/joint-implaus-",combs[1,i], "-", combs[2,i], "-", combs[3,i], ".pdf", sep="")
  pdf(buffer)
  cat("# obs: ", desNames[combs[3,i]], "\n")
  fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE, title.in=paste(mwString, "vs", mwStringCompare, " z: ", desNames[combs[3,i]]))
  dev.off()
}
