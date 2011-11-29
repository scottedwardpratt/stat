## make an analysis using the all the observables we have
## we'll start by combining the metallicity and lum-fn data
## and generating some principle components
##
##
## suppose that we have some R observables
## we can define a multivariate implausibility as I^2 = (E(f(x)) - z )^t V^(-1) (E(f(x)) - z)
## where V is an RxR matrix of Cov(R_i, R_j) with V_obs + V_model added along the diagonal
##
## we can obtain a posterior estimate of  the covariance
## for a pair of functional locations (observations) using the pca decomposition,
## we to rotate the inde Variances of each component back into the real space, this
## estimate will only be as good as the prior "sample" mutual covariance used to create
## the pca decomp, but now we can estimate this at all locations spanned by the design
## 
## need to validate the combined emulator quite a bit before trying to fix the rest, 
## the implausibility becomes very strange when we include contributions from the metallicity
## as well as contributions from the lum.

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

nsamples <- 45  ## we'll hold back 5

## load the model data
##
## first the luminosity data
lumOutputFile  <- "./wave-1/lum_fun_outputs.dat"
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]
##
## now the metallicity data
metOutputFile <- "./wave-1/metallicity_MV_outputs.dat"
modelDataMet <- abs(as.matrix(read.table(metOutputFile)))
nbinsMet <- dim(modelDataMet)[2]

if(nruns != dim(modelDataMet)[1]){
  stop("nruns modelDataMet doesn't match modelDataLum")
}

## redefine nruns
nruns <- nsamples
## now this only works each time we rebuild the model
samp.index <- sample(seq(1,45), size=nsamples)

## load the model data
nbins <-  nbinsMet + nbinsLum 
modelData.big <- cbind(modelDataLum, modelDataMet)
modelData <- modelData.big[samp.index,]

## load the design
designFile <- "./design/design_sorted_wave_1.dat"
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designData.big <- as.matrix(read.table(designFile, col.names=desNames))

designData <- designData.big[samp.index,]

## create the validation set
#holdBack <- list(des=designData.big[-samp.index,], training=modelData.big[-samp.index,])


## load the experimental data
##
## lum
expDataFileLum <- "./lum_fun_observed.dat"
expDataLum <- as.matrix(read.table(expDataFileLum))

expDataFileMet <- "./metallicity_MV_observed.dat"
expDataMet <- abs(as.matrix(read.table(expDataFileMet)))

## the lum data errors are "poisson" ~ sqrt(\lambda), if N_obs was very big we could approximate it as a normal
## distribution with mean \lambda and sd \sqrt{\lambda}
## 
## the met data are 95% confidence bounds, 95% of the data fall within mu +- 2 \sigma
## so sigma = upper-lower / 2 
## where upper = mean + conf, lower = mean - conf
## so sigma = ((mean + conf ) - (mean - conf))  / 2 = conf
expData <- list(obsValue=c(expDataLum[2,], expDataMet[1,]),
                  obsError=c(expDataLum[3,], expDataMet[2,]))
## expData <- list(obsValue=c(expDataLum[2,]),
##                obsError=c(expDataLum[3,]))

rebuild <- 0
buffer <- "functional-data-lum-abs-met.dat"
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

## plot in dimensions 1 and 3 and step in 2
stepData <- fn.emu.steps(fnData, 1, 3, 2)
fn.plot.steps(fnData, stepData,  2)
## you need to run the emulator over steps  to make predictions using fn.emu.steps
## before computing the implausibility
impSteps <- fn.implaus.steps(fnData, stepData)
fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE)

## here we'll plot predictions from the emulator with the dimensions 1,2 and stepping in 3
stepData.escp <- fn.emu.steps(fnData, 1, 2, 3)
impSteps.escp <- fn.implaus.steps(fnData, stepData.escp)
fn.plot.imp.steps(fnData, impSteps.escp, plot.joint=TRUE)



## test how stable the hyperparameters are, create them a bunch of times
## superficially the answer is not very stable but not totally awful
## fnData <- fn.sample.gen(cov.fn=1, reg.order=1)
## fnData <- fn.pca.gen(fnData, cutOff=0.99)

## gen.testThetas <- function(ntest=25){
##   fn.list <- vector("list", ntest)
##   for(i in 1:ntest){
##     fn.list[[i]] <- fn.estimate(fnData)
##   }
##   save(fn.list, file="fndata-thetas-test.dat")
##}
