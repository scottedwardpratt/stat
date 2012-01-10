library("plotrix")
library("calibrate")

##
## ccs, look at the data for the wave-1-long data
## should be about 200 training points, interesting to compare

## plot the luminosity data
## we have the experimentally observed data here
expDataFile <- "lum_fun_observed_2.dat"

## four luminosities have been selected, they are listed in the first row of the file
## the 2nd row contains the values of the cumulative lum fn  (L)
## the 3rd row contains "errors = poisson errors -> sqrt(n)" 
expData <- as.matrix(read.table(expDataFile))
expLumBins <- expData[1,]
expLumValues <- expData[2,]
expLumErrs <- expData[3,]

# the model path for the first wave
modelOutputFile <- "./wave-1-long/lum_fun_outputs_2.dat"
modelData <- as.matrix(read.table(modelOutputFile))
nbins <- dim(modelData)[2]
nruns <- dim(modelData)[1]

# the design file for wave 1 (need to label this better)
designFile <- "./design/design_2_sorted.dat"
desNames <- c("Zr", "Fescp", "Fbary")
# have to correct the design data, it currently includes an "unrun" run # 46 (i think)
designData <- as.matrix(read.table(designFile, col.names=desNames))[1:nruns,]



plotExpVsModel <- function(){
  # plot the experimental data first
  plotCI(expLumBins, expLumValues, uiw=expLumErrs, xlab="luminosity bin", ylab="cumulative luminosity", lwd=2)
  # load the output
  cat("nbins: ", nbins, " \n")
  cat("nruns: ", nruns, " \n")

  for(i in 1:nbins){
    binValue <- rep(expLumBins[i], nruns)    
    points(binValue, modelData[,i], pch=2, col=i)
  }
  plotCI(expLumBins, expLumValues, uiw=expLumErrs, xlab="luminosity bin", ylab="cumulative luminosity", lwd=2, add=TRUE)
  
}

# plot each of the measured observables against all the runs
# aranged by run number so we can get a good feeling for the variation
# we pick the observables by binNumber [1..4] currently
plotExpVsRun <- function(binNumber, pointCol=2){

  buffer <- paste("luminosity bin: ", expLumBins[binNumber], sep="")
  
  plot(modelData[,binNumber], pch=2, col=pointCol,
       xlab="run index", ylab="cumulative lum", main=buffer, type="p")
  plotCI(1, expLumValues[binNumber], uiw=expLumErrs[binNumber], add=TRUE, col=1)
  
  abline(h=expLumValues[binNumber] + expLumErrs[binNumber], lwd=1, lty=2)
  abline(h=expLumValues[binNumber] - expLumErrs[binNumber], lwd=1, lty=2)

}

# plot all the bins using plotExpVsRun
plotAllBinsExpVsRun <- function(){
  par(mfrow=c(2,2)) # this assumes nbins = 2
  for(i in 1:nbins){
    plotExpVsRun(i)
  }
}


plotAllBinsVsAllDes <- function(){
  ndes <- 3
  par(mfrow=c(nbins,ndes), mar=c(4,4,1,1), cex.axis=0.8, cex.lab=0.8, cex.main=0.8, tcl=0.2, mgp=c(2,0.1,0))
#par(mfrow=c(nbins,ndes))
  for(i in 1:nbins){
    for(j in 1:ndes){
      plotBinVsDesign(i, j, 2+i)
    }
  }
}

# plot a bin against a design value
plotBinVsDesign <- function(binNumber, designVariable, pointCol=2){
  buffer <- paste("luminosity bin: ", expLumBins[binNumber], sep="")

  
  if(designVariable == 1){
      plot(designData[,designVariable], modelData[,binNumber], pch=2, col=pointCol,
        xlab=desNames[designVariable], ylab="cumulative luminosity", main=buffer, type="p", las=0, axes=TRUE)
  } else {
    plot(designData[,designVariable], modelData[,binNumber], pch=2, col=pointCol,
        xlab=desNames[designVariable], ylab="", main=buffer, type="p", las=0, axes=TRUE)
  }
      
  abline(h=expLumValues[binNumber], lwd=1, lty=1, col='blue')
  abline(h=expLumValues[binNumber] + expLumErrs[binNumber], lwd=1, lty=2, col='blue')
  abline(h=expLumValues[binNumber] - expLumErrs[binNumber], lwd=1, lty=2, col='blue')

  
}

# make all the pdf plots for checking things
remakePDFS <- function(){
  pdf("./images/lum-exp-vs-model.pdf")
  plotExpVsModel()
  dev.off()

  pdf("./images/lum-all-vs-design.pdf", width=11, height=8.5)
  plotAllBinsVsAllDes()
  dev.off()

  pdf("./images/lum-all-vs-run.pdf")
  plotAllBinsExpVsRun()
  dev.off()

}
