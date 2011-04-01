## computePoints:
## ccs, cec24@phy.duke.edu 2011-March/April (msu)
##
## produces mean and variances for a set of points in the parameter space
## 
## input: a temp file containing the pts produced by the wrapper computePoints.sh
##        also requires the R-object files SampleTemp.dat, PcaTemp.dat and the
##        plain text file theta-table.dat all created by estimateThetas.[sh|R]
## 
## env: rescale=(1|0) rescales and uncenters the output mean and var, or doesn't
##
## output: (stdout) for n input points, 2n lines, of interleaved mean and then variance
## i.e mean(point1) \n var(point1) \n mean(point2) \n var(point2)...
## 
source("gen-samp.R") # load the sample set
source("emu-pca.R") # load the pca fns
# note, need to change this to .so or whatever for deployment
dyn.load("~/local/lib/libRBIND.dylib") #
# moving these files into the project root (or installing them) is probably best
source("~/local/include/libRbind/EmuRbind.R")


## see if we should rescale
cargs <- Sys.getenv(c('rescale', 'rescaleErr', 'inpath'))
rescale <- as.numeric(cargs[1])
rescaleErr <- as.numeric(cargs[2])
inpath <- cargs[3]

## ccs: for debugging
#rescale <- FALSE
#rescaleErr <- TRUE

cat("#rescale = ", rescale, "\n")
cat("#rescaleErr = ", rescaleErr, "\n")
cat("#inpath ", inpath, "\n")

if(rescale == 1){
  cat("#rescaling data with sample mean & var\n")
} else if(rescaleErr == 1){
  cat("#rescaling data with sample mean & prior errors\n")
}


## load up the saved things
## craps out if the files are not there, not ideal
buffer <- paste(inpath, "/SampleTemp.dat", sep="")
load(buffer) # becomes sample
buffer <- paste(inpath, "/PcaTemp.dat", sep="")
load(buffer) # becomes pca.decomp
buffer <- paste(inpath, "/theta-table.dat", sep="")
thetas.est <- read.table(buffer) 

## now we want to read the stdin
#points <- read.table("/dev/stdin") # sadly this doesn't work
points <- read.table("temp")
npoints <- dim(points)[1]

# this is for debug
cat("#output from computePoints.R: mean and then var for each point\n")
cat("#thetas used to generate this in: ./theta-table.dat\n")
cat("#npoints read: ", npoints, "\n")
cat("#ntps: ", pca.decomp$ntps, "\n")

rvalues <- vector("list", npoints)
rvaluesRescaled <- vector("list", npoints)

# compute the curves in t-space at each point
for(i in 1:npoints){
  rvalues[[i]] <- reconCurveAtPoint(points[i,], thetas.est, pca.decomp)
  if(rescale==1){
  # undo the centering and scaling we did in gen.sample
  # with the stored data from sample
    rvaluesRescaled[[i]]$mean <- rvalues[[i]]$mean * sample$sds + sample$means
    rvaluesRescaled[[i]]$var <- rvalues[[i]]$var * sample$sds
    # copy in the unmodified things
    rvaluesRescaled[[i]]$xpos <- rvalues[[i]]$xpos
    rvaluesRescaled[[i]]$tvalues <- rvalues[[i]]$tvalues
  } else if(rescaleErr ==1){
    cat("#using rescale error!\n")
    ## escale with predetermined errors instead of using the sampleerrors
    rvaluesRescaled[[i]]$mean <- (rvalues[[i]]$mean + sample$means) * sample$errs[i,]
    # this is not coming out to the right scale at all
    # \todo fix the rescaled variance
    rvaluesRescaled[[i]]$var <- rvalues[[i]]$var * sample$errs[i,]

    # copy in the unmodified things
    rvaluesRescaled[[i]]$xpos <- rvalues[[i]]$xpos
    rvaluesRescaled[[i]]$tvalues <- rvalues[[i]]$tvalues
  }
}

# now we want to scale the curves back to the values we had before centering
for(i in 1:npoints){
  write.table(rvaluesRescaled[[i]]$mean, stdout(), row.names=FALSE, col.names=FALSE)
  write.table(rvaluesRescaled[[i]]$var, stdout(), row.names=FALSE, col.names=FALSE)
}



