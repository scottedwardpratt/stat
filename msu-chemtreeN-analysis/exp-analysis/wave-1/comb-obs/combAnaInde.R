## combine the metallicity and lum fn information
## but treat each bin as an independent observable (as we've been doing so far)
##
##
## this is silly, we have no idea of the mutual covariance...
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
modelDataMet <- as.matrix(read.table(metOutputFile))
nbinsMet <- dim(modelDataMet)[2]

if(nruns != dim(modelDataMet)[1]){
  stop("nruns modelDataMet doesn't match modelDataLum")
}

nbins <-  nbinsMet + nbinsLum 

modelData <- cbind(modelDataLum, modelDataMet)

## load the design
designFile <- "./design/design_sorted_wave_1.dat"
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designData <- as.matrix(read.table(designFile, col.names=desNames))[1:nruns,]

## load the experimental data
##
## lum
expDataFileLum <- "./lum_fun_observed.dat"
expDataLum <- as.matrix(read.table(expDataFileLum))

expDataFileMet <- "./metallicity_MV_observed.dat"
expDataMet <- as.matrix(read.table(expDataFileMet))

expData <- list(obsValue=c(expDataLum[2,], expDataMet[1,]),
                obsError=c(expDataLum[3,], expDataMet[2,]))

# do we want to rebuild the hyperparameters or not? 1 for yes
rebuild <- 0
buffer <- paste("./thetas-comb-inde.dat")

if(rebuild == 1){
  estimResult <- doCombEstimation(fixNugget=NULL)
  cat("# thetas, scaled design and scaled output saved in: ", buffer, " \n")
  save(estimResult, file=buffer)
} else { # load a previously saved version
  if(file.exists(buffer)){
    cat("# thetas, scaled design and scaled output read from: ", buffer, " \n")
    load(buffer)
  } else {
    estimResult <- doCombEstimation(fixNugget=NULL)
    cat("# thetas, scaled design and scaled output saved in: ", buffer, " \n")
    save(estimResult, file=buffer)
  }    
}



## generate ntest sets of thetas and then dump them to a binary file
## want to see how stable the results are
##
## use the default linear regression model
thetaStabTest <- function(ntest=25, filename="inde-thetas-test-linear.dat"){
  ## we'll use the same list format as the fn data so we can
  ## plot everything the same way...
  fn.list <- vector("list", ntest)
  for(i in 1:ntest){
    temp <- doCombEstimation(fixNugget=NULL)
    fn.list[[i]]$thetas <- temp
    fn.list[[i]]$model.sample$y <- modelData
    fn.list[[i]]$model.sample$des <- designData
    fn.list[[i]]$model.sample$t <- seq(1, nbins)
  }
  save(fn.list, file=filename)
}


thetaStabTest()

#stepPlotDimension(estimResult, 5, 1, 3, 2, plot.var=TRUE)
#stepPlotDimensionImplaus(4,1,3,2, estim.result=estimResult, exp.data=expData)

# plot the first observable against dimensions 1 and 2 while stepping in dimenison 3
emuStepPlots <- function(){
  for(i in 1:nbins){
    buffer <- paste("./images/emu-mean-stepped-obs-",i,"-cov-",estimResult$cov.fn,".pdf", sep="")
    cat("emuStepping: ", buffer, "\n")
    pdf(buffer)
    stepPlotDimension(estimResult,i,1,3,2)
    dev.off()
  }
}

# plot the implausibility against dims 1 and 2
impStepPlots <- function(){
  for(i in 1:nbins){
    buffer <- paste("./images/implaus-stepped-obs-",i,"-cov-", estimResult$cov.fn,".pdf", sep="")
    cat("implausStepping: ", buffer, "\n")
    pdf(buffer)
    stepPlotDimensionImplaus(i,1,3,2, estim.result=estimResult, exp.data=expData)
    dev.off()
  }
}

# plot the implausibility against dims 1 and 2
impStepPlotsComb <- function(){
    buffer <- paste("./images/implaus-stepped-comb-cov-", estimResult$cov.fn,".pdf", sep="")
    cat("implausStepping: ", buffer, "\n")
    pdf(buffer)
    stepPlotDimensionImplausComb(1,3,2, estim.result=estimResult, exp.data=expData)
    dev.off()
}


# do the implausibility sweep for all observables
implausSweep <- function(){
  for(i in 1:nbins){
    buffer <- paste("grid-implaus-lum-bin-", i, ".csv", sep="")
    cat("sweeping: ", buffer, "\n")
    gridImplausSweep(estimResult, expData, i, 1, 2, 3, fname=buffer)
  }
}


# test the estimation of the thetas we have used here
doEstTestLhoodPlot <- function(){
  for(i in 1:nbins){
    buffer <- paste("./images/lhood-compare-bin-", i, "-cov-", estimResult$cov.fn,".pdf", sep="")
    pdf(buffer)
    buffer <- paste("obs bin", i, " ")
    lhoodCompare(estim.result=estimResult, obsIndex=i, titleAdditional=buffer)
    dev.off()
  }
}

doHoldBack <- function(){
  ## this is very useful, it stops things from just totally *dying* when
  ## we get an error
  options(error=utils::recover)
  # we have to pick which observable to index on
  resList <- vector("list", nbins)
  for(i in 1:nbins){
    obsIndex <- i
    resList[[i]] <- HoldBackTest(estimResult$des.scaled, estimResult$train.scaled[,obsIndex], nSets=5)
  }

  resList
}

plotHoldBack <- function(res.list, des.index){
  for(i in 1:nbins){
    if(i == 1){
      plot(estimResult$des.scaled[,des.index], sqrt(res.list[[i]]$deviates), xlab=desNames[des.index], ylab="(E[y] - y) / Sqrt(V[y])", type="p", col=i, pch=2, log="y")
    } else {
      points(estimResult$des.scaled[,des.index], sqrt(res.list[[i]]$deviates), col=i, pch=2+i)
    }
  }
  abline(h=1, col="black", lty=2)
  legend("topright", leg=c("-15.5", "-11.5", "-7.5", "-3.5"), col=c(1,2,3,4), pch=c(2,3,4,5))
}

## compute the ratio of points pooly treated: sqrt(deviation) > sigVal / totalNumber
## useful as a very basic check of the emulator performance across different bins
holdBackRatios <- function(sigVal = 2){
  res.list <- doHoldBack()
  ratios <- rep(NA, nbins)
  for(i in 1:nbins){
    # compare the number of pts which have sd > sigma Val to the total number
    # this is a bit cryptic
    ratios[i] <- length(res.list[[i]]$deviates[sqrt(res.list[[i]]$deviates) > sigVal]) / length(res.list[[i]]$deviates)
  }
  ratios
}
  
holdBackPlots <- function(){
  res.list <- doHoldBack()
  for(i in 1:nparams){
    buffer <- paste("./images/hold-back-", desNames[i],".pdf", sep="")
    pdf(buffer)
    plotHoldBack(res.list, i)
    dev.off()
  }
}
