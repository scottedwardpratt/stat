# lets run the emulator over the luminosity data
# should be able to do this automatically? without changing names between arches?

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

library(MASS)

source("~/local/include/libRbind/emuOverDesign.R") # functions for running the emulator over the design in high dims
source("~/local/include/libRbind/implausOverDesign.R") # functions for computing the implausibility
source("~/local/include/libRbind/testEst.R") # for testing the estimation of thetas

# the model path for the first wave
modelOutputFile <- "./wave-1/lum_fun_outputs.dat"
modelData <- as.matrix(read.table(modelOutputFile))
nbins <- dim(modelData)[2]
nruns <- dim(modelData)[1]

# the design file for wave 1 (need to label this better)
designFile <- "./design/design_sorted_wave_1.dat"
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designDefaults <- c(10, 50, 0.05)
# have to correct the design data, it currently includes an "unrun" run # 46 (i think)
designData <- as.matrix(read.table(designFile, col.names=desNames))[1:nruns,]

expDataFile <- "./lum_fun_observed.dat"
expDataIn <- as.matrix(read.table(expDataFile))
expData <- list( obsValue=expDataIn[2,], obsError=expDataIn[3,])

# do we want to rebuild the hyperparameters or not? 1 for yes
rebuild <- 1
buffer <- paste("./thetas-comb-lum.dat")

if(rebuild == 1){
  estimResult <- doCombEstimation(fixNugget=NULL, cov.fn=1)
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



## repeat all the bits we need to make the figures and tables in the notes
rptForNotes <- function(){

  # re-estimate result
  cat("re estimating thetas\n")
  estimResult <-  doCombEstimation(fixNugget=NULL)
  
  # dump thetas
  thetasFin <- cbind(estimResult$thetas[,1:2], exp(estimResult$thetas[,-(1:2)]))
  write.table(round(thetasFin,3), "thetas-table.dat", row.names=FALSE, col.names=FALSE)

  # do holdbackPlots
  holdBackPlots()
  
  # hold back ratios
  rats <- holdBackRatios()
  write.table(round(rats,3), "hold-back-ratios.dat", row.names=FALSE, col.names=FALSE)
  
  # plot the lhoods
  doEstTestLhoodPlot()

  # plot the emulated means
  emuStepPlots()

  # plot the implausibilities
  impStepPlots()

  #  do the weakly combined implausibility
  impStepPlotsComb()

  # dump a grid swept across the combined imp, for paraview
  gridImplausComb(estim.result=estimResult, exp.data=expData, thresh=1.8, 1, 3, 2)

  # finally rebuild all the pdfs
  source("plotLum.R")
  remakePDFS()
}

# plot the first observable against dimensions 1 and 2 while stepping in dimenison 3
emuStepPlots <- function(){
  for(i in 1:nbins){
    buffer <- paste("./images/emu-mean-stepped-lum-bin-",i,"-cov-",estimResult$cov.fn,".pdf", sep="")
    cat("emuStepping: ", buffer, "\n")
    pdf(buffer)
    stepPlotDimension(estimResult,i,1,3,2)
    dev.off()
  }
}

# plot the implausibility against dims 1 and 2
impStepPlots <- function(){
  for(i in 1:nbins){
    buffer <- paste("./images/implaus-stepped-lum-bin-",i,"-cov-", estimResult$cov.fn,".pdf", sep="")
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







