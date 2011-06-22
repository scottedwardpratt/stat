## testing the problem of the steppiness
nloopPoints <- 50
nparams <- 7
# where we're going to read the emulator data from
inpath <- "../model-data/rhic-2nd-test/"


generatePointList <- function(kVary, nparams, nloopPoints){
  ## fix all the parameters to their central values
  ## then for each param, loop from range_min to range_max
  ## while holding the rest fixed and collect y_0 the first
  ## principle component

  paramRanges <- matrix(0, nrow=nparams, ncol=3)
  ## set the values by hand (read in from ranges.dat in the data folder)
  paramRanges[1,1] <- 0.6 #HYDRO_T0
  paramRanges[1,2] <- 1.1
  paramRanges[2,1] <- 0.6 #GLAUBER_WINBIN_RATIO
  paramRanges[2,2] <- 1.0
  paramRanges[3,1] <- 2.4 #GLAUBER_K_TAU
  paramRanges[3,2] <- 3.3
  paramRanges[4,1] <- 0.5 # HYDRO_INIT_NS
  paramRanges[4,2] <- 2.0
  paramRanges[5,1] <- 0.4 # HYDRO_INIT_FLOW
  paramRanges[5,2] <- 1.2
  paramRanges[6,1] <- 0.03 # HYDRO_SVRATIO
  paramRanges[6,2] <- 0.25
  paramRanges[7,1] <- 0.2 # EQOFST_BUMP_HEIGHT
  paramRanges[7,2] <- 0.8

  for(i in 1:nparams) # set midpoint values (for want of any other idea?)
    paramRanges[i,3] <- 0.5*(paramRanges[i,1]+paramRanges[i,2])

                                        # the list of points we'll evaluate for each loop
  pointList <- matrix(0, nrow=nloopPoints, ncol=nparams)


  for(i in 1:nloopPoints){
    for(k in 1:nparams){ # loop over each parameter
      if(k != kVary){
        pointList[i, k] <- paramRanges[k, 3] # set the default value
      } else {
        pointList[i, k] <- paramRanges[k, 1] + ((paramRanges[k,2]-paramRanges[k,1])/(nloopPoints))*i
      }
    }
  }
  pointList
}          

## now we have to load some shit
source("gen-samp.R")
source("emu-pca.R")
# note, need to change this to .so or whatever for deployment
dyn.load("~/local/lib/libRBIND.so") #
# moving these files into the project root (or installing them) is probably best
source("~/local/include/libRbind/EmuRbind.R")

## load up the saved things
## craps out if the files are not there, not ideal
buffer <- paste(inpath, "/SampleTemp.dat", sep="")
load(buffer) # becomes sample
buffer <- paste(inpath, "/PcaTemp.dat", sep="")
load(buffer) # becomes pca.decomp
buffer <- paste(inpath, "/theta-table.dat", sep="")
thetas.est <- read.table(buffer) 

xvals <- matrix(0, ncol=nparams, nrow=nloopPoints)
yvals <- matrix(0, ncol=nparams, nrow=nloopPoints)
yvalsScaled <- matrix(0, ncol=nparams, nrow=nloopPoints)

reconFlag <- 1


for(k in 1:nparams){
# generate a point list
  pointList <- generatePointList(k, nparams, nloopPoints)
# call the emulator (unprojected)
  if(reconFlag == 0){
    emuValUnProf <- emulateAtListNoProject(pointList, thetas.est, pca.decomp)
    yvals[,k] <- emuValUnProf$mean[,1] # the first y component
    yvalsScaled[,k] <- (yvals[,k] + sample$means[1]) * mean(sample$errs[,1])
  }else {
# lets try projected
    emuValUnProf <- reconCurveAtList(pointList, thetas.est, pca.decomp)
    yvals[,k] <- emuValUnProf$mean[1,] # the first y component
  }
  
  
  xvals[,k] <- emuValUnProf$xpos[,k]
  if(k == 1){
    # plot w/o scale to make it all fit
    plot( yvals[,k]/(max(yvals[,k])-min(yvals[,k])), type="b", col=k)
  } else {
    lines(yvals[,k]/(max(yvals[,k])-min(yvals[,k])) , col=k)
  }
}
