##
## copied from combAna but here we're only going to use the
## lum data, ont the lum+mat included in combAna
## want to try and figure out why adding the mat makes it much
## worse
##
## perhaps need to parcoord the model design to check the "spacing" in some sense

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

#library(MASS)

source("~/local/include/libRbind/emuOverDesign.R") # functions for running the emulator over the design in high dims
source("~/local/include/libRbind/implausOverDesign.R") # functions for computing the implausibility
source("~/local/include/libRbind/testEst.R") # for testing the estimation of thetas

source("fnAnalysis.R")

nsamples <- 45 

## load the model data
##
## first the luminosity data
lumOutputFile  <- "./wave-1/lum_fun_outputs.dat"
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]

##
## now the metallicity data
## metOutputFile <- "./wave-1/metallicity_MV_outputs.dat"
## modelDataMet <- as.matrix(read.table(metOutputFile))
## nbinsMet <- dim(modelDataMet)[2]

if(nruns != dim(modelDataMet)[1]){
  stop("nruns modelDataMet doesn't match modelDataLum")
}

## redefine nruns
nruns <- nsamples
## now this only works each time we rebuild the model
samp.index <- sample(seq(1,45), size=nsamples)

nbins <-  nbinsLum
modelData.big <- modelDataLum
modelData <- modelData.big[samp.index,]

## load the design
designFile <- "./design/design_sorted_wave_1.dat"
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designData.big <- as.matrix(read.table(designFile, col.names=desNames))

designData <- designData.big[samp.index,]



## load the experimental data
##
## lum
expDataFileLum <- "./lum_fun_observed.dat"
expDataLum <- as.matrix(read.table(expDataFileLum))


## the lum data errors are "poisson" ~ sqrt(\lambda), if N_obs was very big we could approximate it as a normal
## distribution with mean \lambda and sd \sqrt{\lambda}
## 
## the met data are 95% confidence bounds, 95% of the data fall within mu +- 2 \sigma
## so sigma = upper-lower / 2 
## where upper = mean + conf, lower = mean - conf
## so sigma = ((mean + conf ) - (mean - conf))  / 2 = conf
expData <- list(obsValue=c(expDataLum[2,]),
                  obsError=c(expDataLum[3,]))


rebuild <- 0
buffer <- "functional-data-lum.dat"
if(rebuild == 1 || file.exists(buffer) == FALSE){
  ##
  ## generate a functional sample from the vars in global sope
  fnData <- fn.sample.gen(cov.fn=1, reg.order=1)
  ## now do the pca decomp
  fnData <- fn.pca.gen(fnData, cutOff=0.99)
  ## estimate the thetas
  fnData <- fn.estimate(fnData)

  save(fnData, file=buffer)
} else {
  load(buffer)
}



## stepData <- fn.emu.steps(fnData, 1, 3, 2, range.Min=0.0, range.Max=1.0)
## #fn.plot.steps(fnData, stepData,  2)

## impSteps <- fn.implaus.steps(fnData, stepData)
## pdf("images-lum/implaus-comb.pdf")
## fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE)
## dev.off()

## plot the pca.decomp eigenvalues
pdf("images-lum/pca-decomp.pdf")
plot(fnData$pca$t, fnData$pca.decomp$ur.h[,1], type="b", ylim=c(-1,1), xlab="obs index", ylab="scaled value")
for(i in 2:nbins){
  lines(fnData$pca$t, fnData$pca.decomp$ur.h[,i], type="b", col=i)
}
title(sub="principle components, including all observables")
legend("topright", paste(1:nbins), lty=rep(1,nbins), col=1:nbins)
dev.off()




## test how stable the hyperparameters are, create them a bunch of times
#fnData <- fn.sample.gen(cov.fn=1, reg.order=1)
#fnData <- fn.pca.gen(fnData, cutOff=0.99)

gen.testThetas <- function(ntest=25){
  fn.list <- vector("list", ntest)
  for(i in 1:ntest){
    fn.list[[i]] <- fn.estimate(fnData)
  }

  save(fn.list, file="fndata-thetas-test.dat")
}
