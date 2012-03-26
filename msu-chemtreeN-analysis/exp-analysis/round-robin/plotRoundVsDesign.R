library("plotrix")
## make plots of the various observables against the design
## include all of the comparison galaxies as a set of exp data

mwlist <- c("mw1", "mw2", "mw3", "mw6")
ncompare <- length(mwlist)

cargs <- Sys.getenv(c('mwAct'))

mwActual <- cargs[1] ## the galaxy that we built the emu for
mwActual.Up <- toupper(mwActual)

                    
cat("mwActual: ", mwActual, "\n")
cat("mwActual.Up: ", mwActual.Up, "\n")


## load the design data
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designFile <- paste("./", mwActual, "/design/design_", mwActual.Up, ".dat", sep="")
designData <- as.matrix(read.table(designFile))
nruns <- dim(designData)[1]

## load the model data
## first the luminosity data
lumOutputFile  <- paste("./", mwActual, "/output/lum_fun_outps_", mwActual.Up, ".dat", sep="")
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]
##
## now the metallicity data
metOutputFile <- paste("./", mwActual, "/output/metallicity_MV_outputs_", mwActual.Up, ".dat", sep="")
modelDataMet <- abs(as.matrix(read.table(metOutputFile)))
nbinsMet <- dim(modelDataMet)[2]

nbins <-  nbinsMet + nbinsLum 
modelData.big <- cbind(modelDataLum, modelDataMet)
modelData <- modelData.big


## now load the exp data for each of the comparison galaxies into its own slice
expDataList <- vector("list", ncompare)
for( index in 1:ncompare){
  expDataFileLum <- paste("./", mwlist[index], "/lum_fun_observed_", toupper(mwlist[index]), ".dat", sep="")
  expDataLum <- as.matrix(read.table(expDataFileLum))

  expDataFileMet <- paste("./", mwlist[index], "/metallicity_MV_observed_", toupper(mwlist[index]), ".dat", sep="")
  expDataMet <- as.matrix(read.table(expDataFileMet))

  expData <- list(obsValue=c(expDataLum[2,], abs(expDataMet[1,])),
                  obsError=c(expDataLum[3,], abs(expDataMet[2,])))

  expDataList[[index]] <- expData
}


## calls all the useful fns in this file
remakePlots <- function(){
  buffer <- paste("./images/", mwActual, "-vs-design.pdf", sep="")
  pdf(buffer, width=10, height=7)
  par(mfrow=c(nbins, nparams), oma=c(4,4,3,4), mar=c(1,4,0,0))
  for(j in 1:nbins){
    if(j == nbins){
      par(mar=c(4,4,0,0))
    }
    plotAgainstDesign(j, color=1)
  }
  par(mfrow=c(1,1))
  title(main=(mwActual), outer=TRUE)
  dev.off()
  
  buffer <- paste("./images/multiplot-", mwActual, ".pdf", sep="")
  pdf(buffer)
  plotAgainstSelf()
  dev.off()

  buffer <- paste("./images/parcoord-", mwActual, ".pdf", sep="")
  pdf(buffer)
  plotPC()
  dev.off()
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




## plot one of the observables against the model
## \todo: needs a legend for the different mw compared
## \todo: needs param labels on the x axis
## \todo: needs obs labels on the y axis
plotAgainstDesign <- function(obsIndex, color){
  cat("# obsIndex ", obsIndex, "\n")
  for(i in 1:nparams){
    plotBinVsDesign(modelData, expData, obsIndex, i, ylab.in=paste(obsIndex), pointCol=color)
  }
}

## plot a bin against a design value
plotBinVsDesign <- function(modelData, expData, binNumber, designVariable, pointCol=2, ylab.in){


  if(designVariable == 1){
      plot(designData[,designVariable], modelData[,binNumber], pch=2, col=pointCol,
        xlab=desNames[designVariable], ylab=ylab.in, type="p", las=0, axes=FALSE)
      axis(2)
  } else {
    plot(designData[,designVariable], modelData[,binNumber], pch=2, col=pointCol,
        xlab=desNames[designVariable], ylab="", type="p", las=0, axes=FALSE)
    axis(2)
  }

  ## put an x axis on the plots at the bottom of the page only
  if(binNumber == nbins){
    axis(1)
  }

  ## plot all the "experimental" comparison mw's
  for(index in 1:ncompare){
    expValue <- expDataList[[index]]$obsValue[binNumber]
    expVar <- expDataList[[index]]$obsError[binNumber]
    abline(h=expValue, lwd=2, lty=1, col=index)
    abline(h=expValue + expVar, lwd=1, lty=3, col=index)
    abline(h=expValue - expVar, lwd=1, lty=3, col=index)
  }
  
}


remakePlots()

