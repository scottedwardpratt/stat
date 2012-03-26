## use the trained emulator (mwActual0) data to generate implausibility grids
## for the observations of mwCompare
##
## ccs, cec24@phy.duke.edu
## 25.01.2012
## copied verbatim from plotImplaus, with different final fns

source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings
initEmu()
source("~/local/include/emu-analysis/fnAnalysis.R")

## 21.02.2012
## assumes that fileNames.R has been sourced
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

expDataLum <- as.matrix(read.table(expDataFileLum))
expDataMet <- as.matrix(read.table(expDataFileMet))

expData <- list(obsValue=c(expDataLum[2,], abs(expDataMet[1,])),
                obsError=c(expDataLum[3,], abs(expDataMet[2,])))

## now we load the fn data
buffer <- paste("functional-data-", nruns, "-test-", mwString, ".dat", sep="")
load(buffer)



## now wodge in the expData, but actually the fn fn.implaus.slice uses the globally available expData
fnData$exp.obs = expData

####
## 
## this is the new bit
##
####

implausGrid <- fn.implaus.grid(fnData, 1, 2, 3, grid.size=32) ## make grid of implausibility in 1,2,3 space

buffer <- paste("./", mwStringCompare , "/implaus-grid-observations-", mwStringCompare, "-training-", mwString, "-dump.csv", sep="")

dump.grid.csv(implausGrid, outname=buffer)

buffer <- paste("./", mwStringCompare , "/implaus-grid-observations-", mwStringCompare, "-training-", mwString, ".dat", sep="")
save(implausGrid, file=buffer)
