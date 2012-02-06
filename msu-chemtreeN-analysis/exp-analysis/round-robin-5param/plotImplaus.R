## use the trained emulator (mwActual0 data to generate implausibility plots
## for the observations of mwCompare

source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings
initEmu()

source("~/local/include/emu-analysis/fnAnalysis.R")
## this defines: {lumOutputFile, metOutputFile, designFile, expDataFileLum, expDataFileMet}
source("./fileNames.R")

## load the model data
##
## first the luminosity data
#lumOutputFile  <- paste("./",mwStringCompare, "/output/lum_fun_outps_", mwString.Up, ".dat", sep="")
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]
##
## now the metallicity data
#metOutputFile <- paste("./", mwStringCompare, "/output/metallicity_MV_outputs_", mwString.Up, ".dat", sep="")
modelDataMet <- abs(as.matrix(read.table(metOutputFile)))
nbinsMet <- dim(modelDataMet)[2]

nbins <-  nbinsMet + nbinsLum 
modelData <- cbind(modelDataLum, modelDataMet)

## load the design
#designFile <- paste("./", mwStringCompare, "/design/design_", mwString.Up, ".dat", sep="")
#desNames <- c("Zr", "Fescp", "Fbary")
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

## generates the joint-implausibility slices
stepData <- fn.emu.steps(fnData, 1, 3, 2)
impSteps <- fn.implaus.steps(fnData, stepData)
buffer <- paste("./", mwStringCompare , "/images/joint-implaus-1-3-2.pdf", sep="")
pdf(buffer)
fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE, title.in=paste(mwString, "vs", mwStringCompare))
dev.off()

## and the other way around
stepData2 <- fn.emu.steps(fnData, 1, 2, 3)
impSteps <- fn.implaus.steps(fnData, stepData2)
buffer <- paste("./", mwStringCompare , "/images/joint-implaus-1-2-3.pdf", sep="")
pdf(buffer)
fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE, title.in=paste(mwString, "vs", mwStringCompare))
dev.off()


