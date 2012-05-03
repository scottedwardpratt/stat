## use the trained emulator (mwActual0) data to generate implausibility grids
## for the observations of mwCompare
##
## ccs, cec24@phy.duke.edu
## 25.01.2012
## copied verbatim from plotImplaus, with different final fns

source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings

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



source("~/local/include/emu-analysis/fnAnalysis.R")


## load the model data
##
## first the luminosity data
lumOutputFile  <- paste("./",mwStringCompare, "/output/lum_fun_outps_", mwString.Up, ".dat", sep="")
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]
##
## now the metallicity data
metOutputFile <- paste("./", mwStringCompare, "/output/metallicity_MV_outputs_", mwString.Up, ".dat", sep="")
modelDataMet <- abs(as.matrix(read.table(metOutputFile)))
nbinsMet <- dim(modelDataMet)[2]

nbins <-  nbinsMet + nbinsLum 
modelData <- cbind(modelDataLum, modelDataMet)



## load the design
designFile <- paste("./", mwStringCompare, "/design/design_", mwString.Up, ".dat", sep="")
desNames <- c("Zr", "Fescp", "Fbary", "sfe", "yde2")
nparams <- length(desNames)
## default point is given as last line, ignore it
designData <- as.matrix(read.table(designFile, col.names=desNames))[1:nruns,]

## now load the observations of the galaxy we want to compare AGAINST
expDataFileLum <- paste("./", mwStringCompare, "/lum_fun_observed_", mwStringCompare.Up, ".dat", sep="")
expDataFileMet <- paste("./", mwStringCompare, "/metallicity_MV_observed_", mwStringCompare.Up, ".dat", sep="")

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
combs <- combn(nparams, 3)
ncombs <- dim(combs)[2]
for(i in 1:ncombs){
  ## this will set fixedValVec to be (0,0,0,0,0) and then the various
  ## projections are made around this point, if you want to project about some
  ## other location you should comment this line
  fixedValVec <- rep(0, nparams)
  ## and uncommment the next one
  ## eg a custom location would be
  ## fixedValVec <- c(0.4, 0.3, 0.1, 0.2, 0.5)
  fixedValVec[combs[1,i]] <- NA
  fixedValVec[combs[2,i]] <- NA
  fixedValVec[combs[3,i]] <- NA
  #stepData <- fn.emu.steps(fnData, combs[1,i], combs[2,i], combs[3,i], fixedValVec)
  #impSteps <- fn.implaus.steps(fnData, stepData)
  #buffer <- paste("./", mwStringCompare, "/images/joint-implaus-",combs[1,i], "-", combs[2,i], "-", combs[3,i], ".pdf", sep="")
  #pdf(buffer)
  #cat("# obs: ", desNames[combs[3,i]], "\n")
  #fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE, title.in=paste(mwString, "vs", mwStringCompare, " z: ", desNames[combs[3,i]]))
  #dev.off()
  implausGrid <- fn.implaus.grid(fnData, combs[1,i], combs[2,i], combs[3,i], fixedValVec=fixedValVec, grid.size=32) ## make grid of implausibility in 1,2,3 space
  buffer <- paste("./", mwStringCompare , "/implaus-grid-observations-",
                  combs[1,i], "-", combs[2,i], "-", combs[3,i], "-mw-",
                  mwStringCompare, "-training-",
                  mwString, "-dump.csv", sep="")
  dump.grid.csv(implausGrid, outname=buffer)
  
  buffer <- paste("./", mwStringCompare , "/implaus-grid-observations-",
                  combs[1,i], "-", combs[2,i], "-", combs[3,i], "-mw-",
                  mwStringCompare, "-training-", mwString,"-inde-",sep="")
  dump.grid.inde.csv(implausGrid, outstub=buffer)

  buffer <- paste("./", mwStringCompare , "/implaus-grid-observations-",
                  combs[1,i], "-", combs[2,i], "-", combs[3,i], "-mw-",
                  mwStringCompare, "-training-", mwString, ".dat", sep="")
  save(implausGrid, file=buffer)

  
}

