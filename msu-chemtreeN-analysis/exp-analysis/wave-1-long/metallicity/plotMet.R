library("plotrix")
library("calibrate")

## plot the meallicity data
## we have the experimentally observed data here
expDataFile <- "metallicity_MV_observed_2.dat"


## the metallicity is modelled as cf1(x) = a11*x + a21
## a11 and a21 are the experimental values we will compare against
## the output of the code
## the errs here are 90% ci's
expData <- as.matrix(read.table(expDataFile))
expBins <- c("a11 - slope", "a21 - intercept")
expValues <- expData[1,]
expErrs <- expData[2,]

# the model path for the first wave
modelOutputFile <- "./wave-1-long/metallicity_MV_outputs_2.dat"
modelData <- as.matrix(read.table(modelOutputFile))
# this is the number of coeffs we'll deal with (bins is an ok name for now)
nbins <- dim(modelData)[2]
# the number of reps we did
nruns <- dim(modelData)[1]

# the design file for wave 1 (need to label this better)
designFile <- "./design/design_2_sorted.dat"
desNames <- c("Zr", "Fescp", "Fbary")
# have to correct the design data, it currently includes an "unrun" run # 46 (i think)
designData <- as.matrix(read.table(designFile, col.names=desNames))[1:nruns,]


# due to the huge difference between c1 and c2 this is not useful at all
plotExpVsModel <- function(){
  # plot the experimental data first
  # 
  plotCI(c(1,2), expValues, uiw=expErrs, xlab="model coeff", ylab="coeff value", lwd=2)
  # load the output
  cat("nbins: ", nbins, " \n")
  cat("nruns: ", nruns, " \n")

  for(i in 1:nbins){
    binValue <- rep(i, nruns)    
    points(binValue, modelData[,i], pch=2, col=i)
  }
  
}

# plot each of the measured observables against all the runs
# aranged by run number so we can get a good feeling for the variation
# we pick the observables by binNumber currently
plotExpVsRun <- function(binNumber, pointCol=2){

  buffer <- paste("coeff: ", expBins[binNumber], sep="")
  
  plot(modelData[,binNumber], pch=3, col=pointCol,
       xlab="run index", ylab="coefficient", main=buffer, type="p")

  plotCI(1, expValues[binNumber], uiw=expErrs[binNumber], add=TRUE, pch=3,col=1)
  abline(h=expValues[binNumber] + expErrs[binNumber], lwd=1, lty=2)
  abline(h=expValues[binNumber] - expErrs[binNumber], lwd=1, lty=2)

}

# plot all the bins using plotExpVsRun
plotAllBinsExpVsRun <- function(){
  par(mfrow=c(1,2)) # this assumes nbins = 2
  for(i in 1:nbins){
    plotExpVsRun(i)
  }
}


plotAllBinsVsAllDes <- function(){
  ndes <- 3
  par(mfrow=c(nbins,ndes))
  for(i in 1:nbins){
    for(j in 1:ndes){
      plotBinVsDesign(i, j, i)
    }
  }
}

# plot a bin against a design value
plotBinVsDesign <- function(binNumber, designVariable, pointCol=2){

  buffer <- paste("coeff: ", expBins[binNumber], sep="")
  
  plot(designData[,designVariable], modelData[,binNumber], pch=3, col=pointCol,
       xlab=desNames[designVariable],ylab="coefficient", main=buffer, type="p")

   
  abline(h=expValues[binNumber], lwd=1, lty=1, col='blue')
  abline(h=expValues[binNumber] + expErrs[binNumber], lwd=1, lty=2, col='blue')
  abline(h=expValues[binNumber] - expErrs[binNumber], lwd=1, lty=2, col='blue')

  
}

remakePDFS <- function(){
  pdf("images/met-all-vs-run.pdf", width=11, height=7)
  plotAllBinsExpVsRun()
  dev.off()

  pdf("images/met-all-vs-design.pdf")
  plotAllBinsVsAllDes()
  dev.off()

}
