
## rescales the given set of results from computePoints and dumps them to
## stdout
## with std rescale we don't fold back in the original errors (kept in sample)
## with rescaleErr we do fold this error back in (which is what we ought to use)
##
## but this is confusing, we are generating some nEmuPts points and we're trying to rescale
## them with values generated for nDesignPoints
## need to be careful here...
##
rescaleAndOuputRvalues <- function(rvalues, fn.data, npoints){
  trainCenter.vec <- attr(fn.data$model.sample$y, "scaled:center")
  trainScale.vec <- attr(fn.data$model.sample$y, "scaled:scale")
  nbins <- dim(fn.data$model.sample$y)[2]

  rvalues.rescaled <- list(mean=matrix(0, nrow=npoints, ncol=nbins),
                           var=matrix(0, nrow=npoints, ncol=nbins))

  rvalues.rescaled$mean <- rvalues$mean*trainScale.vec + trainCenter.vec
  rvalues.rescaled$var <- rvalues$var*(trainScale.vec)^2

  for(i in 1:npoints){
    ## this just to make the output easier to look at
    cat("\n")
    write.table(t(rvalues.rescaled$mean[,i]), stdout(), row.names=FALSE, col.names=FALSE)
    cat("\n")
    write.table(t(rvalues.rescaled$var[,i]), stdout(), row.names=FALSE, col.names=FALSE)

    if(any(is.na(rvalues.rescaled$mean[,i])==TRUE)||any(is.na(rvalues.rescaled$var[,i])==TRUE) ){
      stop("NA values in emulator output.")
    }
    
  }
}

## just a plainish dump to stdout
outputUnProjected <- function(rvalues, npoints){
  cat("########################################################################\n")
  cat("# unreconstructed data follows. Mean and variance alternate lines.     #\n")
  cat("########################################################################\n")
  for(i in 1:npoints){
    write.table(t(rvalues$mean[i,]), stdout(), row.names=FALSE, col.names=FALSE)
    cat("\n")
    write.table(t(rvalues$var[i,]), stdout(), row.names=FALSE, col.names=FALSE)
    cat("\n")
    if(any(is.na(rvalues$mean[i,])==TRUE)||any(is.na(rvalues$var[i,])==TRUE) ){
      stop("NA values in emulator output.")
    }
  }
}

