########################################################################
# ccs: cec24@phy.duke.edu  feb-2011
# england, east-lansing
########################################################################
#
# purpose:
# load up the design and put together a sample-list in the style of the 
# lists used in emu-pca
#
# love is all from what i've heard,
# but my hearts learnt to kill,
# oh mine has learnt to kill


## generate a list which is compatible with the previously
## designed emulator / pca methods
gen.samp <- function(
                     fname = "", # savey?
                     desPath = "", # alt path for design
                     sampPath = "", # alt path for samples
                     sampErrPath = "" #path for errors
                     ){

  if(nchar(desPath) > 0){
    des <- load.design(desPath)
  } else{
    des <- load.design()           # load the default
  }

  if(nchar(sampPath) > 0){
    samp <- load.sample(samPath, sampErrPath)
  } else{
    samp <- load.sample()               # load the default,
                                        #it will be centered
  }

  tvec <- seq(1, samp$ntps) 

  # we could do some error checking here to make sure that nreps and
  # ntps etc all match up
  if(des$nreps != samp$nreps)
    stop("sample and design have differnet number of reps!")

  # the final sample
  # note that the existing pca apparaturs (emu-pca.R)
  # expects the matrices to be size x reps
  # whereas the above defns return matrices which are reps x size
  # hence the transpose
  sample <- list(t=tvec, y=t(samp$y), des=t(des$des))

  if(nchar(fname)>0) {
    save(sample, file=fname);
    cat(paste("Object `", "sample", "' saved in file: `",
                fname,"'.\n", sep=""));
  }

  invisible(sample)
}
                     




## load the design from the parameters folder and assemble it into
## a useful list
load.design <-  function(path="../parameters/stats.param-combined-s"){
  des <- as.matrix(read.table(path, header=FALSE))
  # figure out how many variables are in the design
  # the first column is always the run id so we don't care about that
  ndesVars <- (dim(des)[2] - 1)
  nReps <- dim(des)[1] # the number of repeats
  # and this is what we return
  # nb1: we deliberately cut out the id column
  # nb2: when cutting out the column the final col is nvars+1
  design <- list(nvars=ndesVars, nreps=nReps, des=des[,2:(ndesVars+1)])

  invisible(design)
}


load.sample.old <- function(path="../analysis/combined/analysis.b3.29-combined-s", center=TRUE){
  samp <- read.table(path)
  nReps <- dim(samp)[1]
  nTps <-  (dim(samp)[2] -1)

  sampCut <- samp[,2:(nTps+1)]
  if(center){
    sampCutScaled <- scale(sampCut)
  }
  sample <- list(nreps=nReps, ntps=nTps, y=sampCutScaled)
}

## read combined samples and assemble into a
## useful list
## this will also center the sample data by default
#
# if impacts is true then the second column of the read data is taken to be strings giving
# the impact paramter for each row, this is a PARAMETER and not a result and so we do not include
# it in the sample y's
#
# the sample files need to be sorted by ID already as i'm not smart enough to do this in R
load.sample <- function(path="../analysis/spectra-combined-s.dat",
                        errpath="../analysis/spectra-combined-errs-s.dat",
                        center=TRUE, centerWithErrs=FALSE){

  samp <- read.table(path, header=FALSE)
  errSamp <- read.table(errpath, header=FALSE)
  # now each row is one design rep,
  nReps <- dim(samp)[1]
  #
  # the width of each row is the number of points in the "T-space" i.e
  # the functional space which we're going to try and reduce using
  # the pca
  # 
  # again the first column is the runId so we don't care about this too
  # much, since from here on in we can be sure that R won't mess with
  # the row ordering unless we tell it to.
  #
  nTps <- (dim(samp)[2] - 1) 

  sampCut <- samp[,2:(nTps+1)]
  errSampCut <- errSamp[,2:(nTps+1)]
  
  if(center) {
    ## from the docs
    ## 'scale' is generic function whose default method centers and/or
    ## scales the columns of a numeric matrix
    ## scale(x, center=TRUE, scale=TRUE)
    ##
    ## we need to keep the column means and sds so that we can reconstruct
    ## stuff later. Yscaled = (Y - <YcolMean>) / (YcolSd)
    Means <- rep(0, nTps)
    Sds <- rep(0, nTps)
    errSampCutScaled <- matrix(nrow=nReps, ncol=nTps)
    for(i in 1:nTps){
      Means[i] <- mean(sampCut[,i])
      Sds[i] <- sd(sampCut[,i])
      # scale the errors by the sd's so that they remain the same magnitude when we're using the scaled Y's
      errSampCutScaled[,i] <- errSampCut[,i] / Sds[i]
    }
 
    sampScaled <- scale(sampCut)
  } else if(centerWithErrs) {
    ## here the sample Yp =  (Y / err_y)  - < Y / err_y > 
    Means <- rep(0, nTps)
    ## divide all the data by the given errors and then center (but don't scale)
    Means <- mean(sampCut/errSamp[,2:(nTps+1)])
    sampScaled <- scale(sampCut / errSamp[,2:(nTps+1)], center=TRUE, scale=FALSE)
  } else {
    sampScaled <- sampCut
  }

  if(!center && !centerWithErrs){
    sample <- list(nreps = nReps, ntps = nTps, y=sampScaled, errs=errSampCut )
  } else if (center) {
    # we can reconstruct any of the y-values by doing
    # yrecon[index] = yscaled[,index] * sds[index] + means[index]
    sample <- list(nreps = nReps, ntps = nTps, y=sampScaled, means=Means, sds=Sds, errs=errSampCutScaled)
  } else if(centerWithErrs){
    sample <- list(nreps = nReps, ntps = nTps, y=sampScaled, means=Means, errs=errSampCut)
  }
  invisible(sample)
}



## test this all out
#test <- gen.samp(fname="sample-spec-test.dat")

  
  

