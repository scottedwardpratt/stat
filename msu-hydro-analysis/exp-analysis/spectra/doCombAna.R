## analyse the spectrum data
source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings
# now load the emulator lib
initEmu()
## load the functional analysis routines
source("~/local/include/emu-analysis/fnAnalysis.R")

## load the shared information for this run
source("./setup.R")


## load the design data
nparams <- 6
designData <- as.matrix(read.table(designFile))[,2:(nparams+1)]
desNames <- designFieldNames

modelData <- as.matrix(read.table(modelDataFile, header=TRUE))
colNames <- dimnames(modelData)[[2]]
nbins <- dim(modelData)[2]
nruns <- dim(modelData)[1]

expData.in <- t(as.matrix(read.table(defaultFile, header=TRUE)))
expData <- list(obsValue=c(as.numeric(expData.in[,1])), obsError=(as.numeric(expData.in[,2])))

rebuild <- 0
buffer <- paste("functional-data-", analysisBaseName, "thetas.dat", sep="")
if(rebuild == 1 || file.exists(buffer) == FALSE){
  ##
  ## generate a functional sample from the vars in global sope
  fnData <- fn.sample.gen(cov.fn=1, reg.order=0)
  ## now do the pca decomp
  fnData <- fn.pca.gen(fnData, cutOff=0.99)
  ## estimate the thetas
  fnData <- fn.estimate(fnData)
  
  save(fnData, file=buffer)
} else {
  load(buffer)
}


## systematically generate all projections of the joint implaus
combs <- combn(nparams, 3)
ncombs <- dim(combs)[2]
for(i in 1:ncombs){
  fixedValVec <- rep(0, nparams)
  fixedValVec[combs[1,i]] <- NA
  fixedValVec[combs[2,i]] <- NA
  fixedValVec[combs[3,i]] <- NA
  stepData <- fn.emu.steps(fnData, combs[1,i], combs[2,i], combs[3,i], fixedValVec)
  impSteps <- fn.implaus.steps(fnData, stepData)
  pdf(paste("./images/implaus-",combs[1,i], "-", combs[2,i], "-", combs[3,i], ".pdf", sep=""))
  fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE)
  dev.off()
}


