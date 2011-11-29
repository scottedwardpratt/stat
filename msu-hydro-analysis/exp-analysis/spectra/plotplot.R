library("plotrix")

## load the shared stuff for this analysis 
source("setup.R")


## load the design data
nparams <- 6 ## now we only have 6 params
designData <- as.matrix(read.table(designFile))[,2:(nparams+1)]
desNames <- designFieldNames
nruns <- dim(designData)[1]


modelData <- as.matrix(read.table(modelDataFile, header=TRUE))
colNames <- dimnames(modelData)[[2]]
nbins <- dim(modelData)[2]

expData <- as.matrix(read.table(defaultFile, header=TRUE))


remakePlots <- function(){
  for(j in 1:nbins){
    buffer <- paste("./images/", analysisBaseName, colNames[j], "-vs-design.pdf", sep="")
    pdf(buffer)
    plotAgainstDesign(j)
    dev.off()
  }

  buffer <- paste("./images/multiplot-", analysisBaseName, ".pdf", sep="")
  pdf(buffer)
  plotAgainstSelf()
  dev.off()

  buffer <- paste("./images/parcoord-", analysisBaseName, ".pdf", sep="")
  pdf(buffer)
  plotPC()
  dev.off()
}

## plot one of the observables against the model
## for a given impact
plotAgainstDesign <- function(obsIndex){
  cat("# obsIndex ", obsIndex, "\n")
  par(mfrow=c(2,4), oma=c(1,1,3,0), mar=c(4,4,3,1))
  for(i in 1:nparams){
    plotBinVsDesign(modelData, expData, obsIndex, i, ylab.in=colNames[obsIndex])
  }
  par(mfrow=c(1,1))
  title(paste(colNames[obsIndex], " cclass=", cclass), outer=TRUE)
}

plotAgainstSelf <- function(){
  pairs(modelData)
}  


plotPC <- function(scale=FALSE){
  library("MASS")
  if(scale==FALSE){
    parcoord(modelData, var.label=TRUE)
  } else {
    parcoord(scale(modelData), var.label=TRUE)
  }
}
  

# plot a bin against a design value
plotBinVsDesign <- function(modelData, expData, binNumber, designVariable, pointCol=2, ylab.in){

  if(designVariable == 1){
      plot(designData[,designVariable], modelData[,binNumber], pch=2, col=pointCol,
        xlab=desNames[designVariable], ylab=ylab.in, type="p", las=0, axes=TRUE)
  } else {
    plot(designData[,designVariable], modelData[,binNumber], pch=2, col=pointCol,
        xlab=desNames[designVariable], ylab="", type="p", las=0, axes=TRUE)
  }

  expValue <- as.numeric(expData[1,binNumber])
  expVar <- as.numeric(expData[2,binNumber])

  abline(h=expValue, lwd=1, lty=1, col='blue')
  abline(h=expValue + expVar, lwd=1, lty=2, col='blue')
  abline(h=expValue - expVar, lwd=1, lty=2, col='blue')

  
}



