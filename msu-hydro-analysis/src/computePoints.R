## computePoints:
## ccs, cec24@phy.duke.edu 2011-March/April/November (msu)
## 1) updated to use reconCurveAtList which is much faster
## than reconCurveAtPoint.
##
## produces mean and variances for a set of points in the parameter space
## 
## input: a temp file containing the pts produced by the wrapper computePoints.sh
##        also requires the R-object files SampleTemp.dat, PcaTemp.dat and the
##        plain text file theta-table.dat all created by estimateThetas.[sh|R]
## 
##
## output: (stdout) for n input points, 2n lines, of interleaved mean and then variance
## i.e mean(point1) \n var(point1) \n mean(point2) \n var(point2)...
##

emulateAtListNoProject <- function(pointList, fnData){
  ## emulate the code at the points given by pointList
  nr <- dim(fnData$pca.decomp$zmat)[1]
  nparams <- dim(fnData$pca.decomp$des)[1]
  nmodelpts <- dim(fnData$pca.decomp$des)[2]
  # points = nemupts x nparams points
  nemupts <- dim(pointList)[1]
  
  thetas <- fnData$thetas.est
  nthetas <- dim(thetas)[2]

  emu.results <- list(mean=matrix(0, ncol=nr, nrow=nemupts),
                      var=matrix(0, ncol=nr, nrow=nemupts))
  
  for(i in 1:nr){
    model <- list(xmodel=as.matrix(t(fnData$pca.decomp$des)), training=t(fnData$pca.decomp$zmat)[,i])

    ## 
    temp <- callEmulateAtList(model, as.numeric(thetas[i,]), as.numeric(pointList),
                              nemupts, nmodelpts, nparams=nparams, nthetas=nthetas,
                              cov.fn=fnData$cov.fn, reg.order=fnData$reg.order)

    emu.results$mean[,i] <- temp$mean
    emu.results$var[,i] <- temp$var

  }
  emu.results
}




cargs <- Sys.getenv(c('inpath','fndata', 'reconstruct','EMU_DIRECTORY'))

inpath <- cargs[1]
fndata.name <- cargs[2]
reconstruct <- as.numeric(cargs[3])
directory <- cargs[4]

## load up the saved things
load(fndata.name)


cat("#inpath: ", inpath, "\n")
cat("#fndata: ", fndata.name, "\n")
cat("#directory: ", directory, "\n")

cat("########################################################################\n")
cat("# reconstruct: ", reconstruct, "                                                     #\n")

if(reconstruct == 1){
  cat("# projecting data back into functional space\n")
} else {
  cat("# data is *not* projected                                              #\n");
  cat("# nr: ", fnData$pca.decomp$nr, "                                                              #\n")
  cat("########################################################################\n")
}

## moving these files into the project root (or installing them) is probably best
#source("computePointsOutput.R")
buffer <- paste(directory, "/computePointsOutput.R", sep="")
source(buffer)

source("~/local/include/libRbind/EmuRbind.R")
    ## load the functional analysis routines
source("~/local/include/emu-analysis/fnAnalysis.R")

## now load the emu library
#initEmu()

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


## now we want to read the stdin
points <- as.matrix(read.table("temp"))
npoints <- dim(points)[1]
numpoints <-dim(points)[2]

# this is for debug
cat("# output from computePoints.R: mean and then var for each point\n")
cat("# npoints read: ", npoints, "\n")
cat("# ntps: ", fnData$pca.decomp$ntps, "\n")
cat("########################################################################\n")

rvalues <- vector("list", npoints)

if(reconstruct == 1){
# compute the curves in t-space at each point
  rvalues <- reconCurveAtList(points, fnData$thetas.est, fnData$pca.decomp,
                              cov.fn=fnData$cov.fn, reg.order=fnData$reg.order)
} else {
  rvalues <- emulateAtListNoProject(points, fnData)
}

if(reconstruct == 1) { # if we have done the projection we need to rescale the data
  rescaleAndOuputRvalues(rvalues, fnData, npoints);
} else{
  outputUnProjected(rvalues, npoints)
}



