## computePoints:
## ccs, cec24@phy.duke.edu 2011-March/April (msu)
## 1) updated to use reconCurveAtList which is much faster
## than reconCurveAtPoint.
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

## see if we should rescale
cargs <- Sys.getenv(c('rescale', 'rescaleErr', 'inpath', 'reconstruct','EMU_DIRECTORY'))
rescale <- as.numeric(cargs[1])
rescaleErr <- as.numeric(cargs[2])

inpath <- cargs[3]
reconstruct <- as.numeric(cargs[4])
directory <- cargs[5]

cat("#inpath: ", inpath, "\n")
cat("#directory: ", directory, "\n")

cat("########################################################################\n")
cat("# reconstruct: ", reconstruct, "                                                     #\n")

if(reconstruct == 1){
  cat("# projecting data back into functional space\n")
} else {
  cat("# data is *not* projected                                              #\n");
  cat("########################################################################\n")
}

source(paste(directory, "/gen-samp.R", sep = "")) # load the sample set
source(paste(directory, "/emu-pca.R", sep = "")) # load the pca fns


# note, need to change this to .so or whatever for deployment
dyn.load("~/local/lib/libRBIND.so") #
# moving these files into the project root (or installing them) is probably best
source("~/local/include/libRbind/EmuRbind.R")
source("computePointsOutput.R")
#source("~/local/include/libRbind/multivar.R")




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
points <- as.matrix(read.table("temp"))
npoints <- dim(points)[1]
numpoints <-dim(points)[2]

# this is for debug
cat("#output from computePoints.R: mean and then var for each point\n")
cat("#thetas used to generate this in: ./theta-table.dat\n")
cat("#npoints read: ", npoints, "\n")
cat("#ntps: ", pca.decomp$ntps, "\n")

rvalues <- vector("list", npoints)
rvaluesRescaled <- vector("list", npoints)

if(reconstruct == 1){
# compute the curves in t-space at each point
# this version takes < 1 second for 330 points,
# much more efficient.
  rvalues <- reconCurveAtList(points, thetas.est, pca.decomp)
# REFERENCE ONLY
# this was the old method, a lot of setup and teardown had to be done
# for each point, making it very slow
## print(points)
## for(i in 1:336){ # 17 secondsp
## this isn't actually practical, it's just to test how much slower a single point would be
##   rvaluesRescaled <- reconCurveAtPoint(points[1,], thetas.est, pca.decomp)
## }
} else {
  rvalues <- emulateAtListNoProject(points, thetas.est, pca.decomp)
}

if(reconstruct == 1) { # if we have done the projection we need to rescale the data
  rescaleAndOuputRvalues(rvalues, rvaluesRescaled, sample, rescale, rescaleErr, npoints);
} else{
  outputUnProjected(rvalues, npoints)
}





