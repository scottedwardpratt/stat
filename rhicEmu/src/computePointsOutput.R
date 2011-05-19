
## rescales the given set of results from computePoints and dumps them to
## stdout
## with std rescale we don't fold back in the original errors (kept in sample)
## with rescaleErr we do fold this error back in (which is what we ought to use)
##
rescaleAndOuputRvalues <- function(rvalues, rvaluesRescaled, sample, rescale=0, rescaleErr=1, npoints){
  for(i in 1:npoints){
    if(rescale==1){
                                        # undo the centering and scaling we did in gen.sample
                                        # with the stored data from sample
      rvaluesRescaled[[i]]$mean <- rvalues$mean[,i] * sample$sds + sample$means
      rvaluesRescaled[[i]]$var <- rvalues$var[,i] * sample$sds
                                        # copy in the unmodified things
      rvaluesRescaled[[i]]$xpos <- rvalues$xpos[i,]
      rvaluesRescaled[[i]]$tvalues <- rvalues$tvalues
    } else if(rescaleErr ==1){
      cat("#using rescale error!\n")
      ## print(rvalues$mean[,i])
      ## print(rvalues$xpos[i,])
      ## escale with predetermined errors instead of using the sampleerrors
      rvaluesRescaled[[i]]$mean <- (rvalues$mean[,i] + sample$means) * sample$errs[i,]
                                        # this is not coming out to the right scale at all
                                        # \todo fix the rescaled variance
                                        # errors combine in quadrature, dumbass
      rvaluesRescaled[[i]]$var <- sqrt((rvalues$var[,i]*rvalues$var[,i]) + (sample$errs[i,]*sample$errs[i,]))

                                        # copy in the unmodified things
      rvaluesRescaled[[i]]$xpos <- rvalues$xpos[i,]
      rvaluesRescaled[[i]]$tvalues <- rvalues$tvalues
    }
  }

                                        # now we want to scale the curves back to the values we had before centering
  for(i in 1:npoints){
                                        # this just to make the output easier to look at
    cat("\n")
    write.table(rvaluesRescaled[[i]]$mean, stdout(), row.names=FALSE, col.names=FALSE)
    cat("\n")
    write.table(rvaluesRescaled[[i]]$var, stdout(), row.names=FALSE, col.names=FALSE)
    if(any(is.na(rvaluesRescaled[[i]]$mean)==TRUE)||any(is.na(rvaluesRescaled[[i]]$var)==TRUE) ){
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

