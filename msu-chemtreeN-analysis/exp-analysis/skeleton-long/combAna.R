## make an analysis using the all the observables we have
## we'll start by combining the metallicity and lum-fn data
## and generating some principle components
##
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

##
## the 200 point data set ends up with particularly small variances 
## compared to the 45 point set.
##
## should perhaps try doing a 100 and a 150 point set too?
## 

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

nsamples <- 200  ## we'll hold back 10

## load the model data
##
## first the luminosity data
lumOutputFile  <- "./wave-1-long/lum_fun_outputs_2.dat"
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]
##
## now the metallicity data
metOutputFile <- "./wave-1-long/metallicity_MV_outputs_2.dat"
modelDataMet <- abs(as.matrix(read.table(metOutputFile)))
nbinsMet <- dim(modelDataMet)[2]

if(nruns != dim(modelDataMet)[1]){
  stop("nruns modelDataMet doesn't match modelDataLum")
}

## redefine nruns
nruns <- nsamples


nbins <-  nbinsMet + nbinsLum 
#nbins <- nbinsLum
modelData.big <- cbind(modelDataLum, modelDataMet)
#modelData.big <- modelDataLum

modelData <- modelData.big

#samp.index <- sample(seq(1,200), size=nsamples)
#modelData <- modelData.big[samp.index,]


## load the design
designFile <- "./design/design_2_sorted.dat"
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designData.big <- as.matrix(read.table(designFile, col.names=desNames))

designData <- designData.big

## create the validation set
##holdBack <- list(des=designData.big[-samp.index,], training=modelData[-samp.index,])

## load the experimental data
##
## lum
expDataFileLum <- "./lum_fun_observed_2.dat"
expDataLum <- as.matrix(read.table(expDataFileLum))

expDataFileMet <- "./metallicity_MV_observed_2.dat"
expDataMet <- as.matrix(read.table(expDataFileMet))

## the lum data errors are "poisson" ~ sqrt(\lambda), if N_obs was very big we could approximate it as a normal
## distribution with mean \lambda and sd \sqrt{\lambda}
## 
## the met data are 95% confidence bounds, 95% of the data fall within mu +- 2 \sigma
## so sigma = upper-lower / 2 
## where upper = mean + conf, lower = mean - conf
## so sigma = ((mean + conf ) - (mean - conf))  / 2 = conf
expData <- list(obsValue=c(expDataLum[2,], abs(expDataMet[1,])),
                  obsError=c(expDataLum[3,], abs(expDataMet[2,])))


rebuild <- 1
buffer <- paste("functional-data-lum-", nsamples, "-test.dat", sep="")
if(rebuild == 1 || file.exists(buffer) == FALSE){
  ##
  ## generate a functional sample from the vars in global sope
  ## regression with a constant term seems to work better in general than
  ## a first order term
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
#fn.plot.steps(fnData, stepData,  2)
## you need to run the emulator over steps  to make predictions using fn.emu.steps
## before computing the implausibility
impSteps <- fn.implaus.steps(fnData, stepData)
#fn.plot.imp.steps(fnData, impSteps, 5, plot.joint=FALSE)
pdf("images/joint-implaus-1-3-2-repeat.pdf")
fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE)
dev.off()


stepData <- fn.emu.steps(fnData, 1, 2, 3)
impSteps <- fn.implaus.steps(fnData, stepData)
pdf("images/joint-implaus-1-2-3-repeat.pdf")
fn.plot.imp.steps(fnData, impSteps, plot.joint=TRUE)
dev.off()




## fnData <- fn.sample.gen(cov.fn=1, reg.order=0)
## fnData <- fn.pca.gen(fnData, cutOff=0.99)
## gen.testThetas <- function(ntest=25){
##   fn.list <- vector("list", ntest)
##   for(i in 1:ntest){
##     fn.list[[i]] <- fn.estimate(fnData)
##   }

##   save(fn.list, file="fndata-thetas-test-200.dat")
## }



plotEmu.steps <- function(){
  for(i in 1:nbins){
    buffer <- paste("./images/", nsamples, "/emu-mean-obs-",i,".pdf", sep="")
    pdf(buffer)
    fn.plot.steps(fnData, stepData, i)
    dev.off()
    buffer <- paste("./images/", nsamples, "/emu-var-obs-",i,".pdf", sep="")
    pdf(buffer)
    fn.plot.steps(fnData, stepData, i, plot.var=TRUE)
    dev.off()
  }
}
    
#plotEmu.steps()
