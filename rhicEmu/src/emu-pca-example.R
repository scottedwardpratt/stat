#####################################################
# Test function: damped sine wave
#####################################################
yM <- function(t, x1, x2)
  { return(exp(-outer(t, x1))*sin(outer(t,12*x2))); }

#####################################################
# Function to generate example samples
#####################################################
gen.samp.artificial <- function(     
      nreps  = 15,             # |I_S| in notes: no. of samples
      ntps   = 100,            # q     in notes: no. of "t" loc'ns
      x1.lim = c(0.20, 0.25),  # Range for vble x1
      x2.lim = c(0.20, 0.25),  # Range for vble x2
      t.lim  = c(0,5),         # Range for "t"
      fname = ""               # Wanna save it?
      ) {

  xlim  <- cbind(x1=x1.lim, x2=x2.lim);
  des   <- t(maximinLHS(nreps, 2) %*% diag(apply(xlim,2,diff)) +
           xlim[rep(1,nreps),]);
  tvec  <- seq(t.lim[1], t.lim[2], , ntps);
  samp  <- yM(tvec, des[1,], des[2,]);
  ## plot(tvec, samp[,1], type="l", col=1,
  ##      xlab="t", ylab=expression(Y^M ~ (t)));
  ## for(i in 2:nreps){
  ##   lines(tvec, samp[,i], col=i);
  ## }
  sample <- list(t=tvec, y=samp, des=des);
  ## x11()
  if(nchar(fname)>0) {
    save(sample, file=fname);
    cat(paste("Object `", "sample", "' saved in file: `",
                fname,"'.\n", sep=""));
  }
  invisible(sample);
}



doZEmulator <- function(design, thetas, pca.decomp){
  ##
  ## compute the mean and variance of a GP with the given hyperparams thetas in the
  ## pca space.
  ##
  ## design: the locations in the parameter space  d_1 = (x1,x2,theta1,theta2,etc) that we will
  ## evaluate the emulated mean and variance at. The min/max range of these values will 
  ## set the coverage of the reconstructed emulated curves, relative to the original set of 
  ## design parameters used to generate their original samples
  ##
  ## thetas: the estimated hyperparameters for a GP trained from the samples of the data
  ## as transformed into the pca space and the original design.
  ##
  ## pca.decomp: contains all the info (generated from emu.pca)
  ##
 
 #thetas.est <- multidim(t(pca.decomp$des), pca.decomp$nreps, t(pca.decomp$zmat), dim(pca.decomp$zmat)[1])
  nr <- dim(pca.decomp$zmat)[1]
  nemupts <- dim(design)[1]
  nparams <- dim(pca.decomp$des)[1]
  nmodelpts <- dim(pca.decomp$des)[2]
  result <- list(des=design, emuMean=matrix(0, nrow=nemupts, ncol=nr), emuVar=matrix(0, nrow=nemupts, ncol=nr))
  #browser()
  for(i in 1:nr){
    model <- list(xmodel=as.matrix(t(pca.decomp$des)), training=t(pca.decomp$zmat)[,i])
    for(j in 1:nemupts){
      # call the emulator at the jth point in the design space
      temp <- callEmulateAtPoint(model, as.numeric(thetas[i,]), as.numeric(design[j,]), nmodelpts, nparams=nparams, nthetas=(nparams+2))
      if(is.nan(temp$mean) | abs(temp$mean) > 1E3){
        result$emuMean[j,i] <- 0
      } else {
        result$emuMean[j,i] <- temp$mean
      }
      result$emuVar[j,i] <- temp$var
    }
  }
  invisible(result)
}

#write.table("zmat.dat", zemu.points)


  

checkZEmuResult <- function(pca.decomp, zemu.result, index){
  ## this requires the scatterplot3d lib
  ## just plot the original data for the first pca cpt in the first two dims (sorry math i can
  ## only just about vis a 2d fn) and then try and draw the emulated surface over it
  library("scatterplot3d")

  ## plot the first 2 cpts of the  original design points as xy coords and the #
  ## pca response in this projectino as the z height
  s3d<-scatterplot3d(pca.decomp$des[1,], pca.decomp$des[2,], pca.decomp$zmat[index,], grid=TRUE, highlight.3d=TRUE, type="h", pch=16)
  ## now do the same for the emulated curves
  s3d$points3d(zemu.result$des[,1], zemu.result$des[,2], zemu.result$emuMean[,index], col="blue", type="h", pch=16)
  title(main="Check the Performance of the Emulator", sub="the blue surface should approximately match up with the black points")
}


reconstructZEmuResult <- function(pca.decomp, zemu.result, emuDesign){
  ## takes the emulated mean and variance in the pca space and reconstructs it
  ## back into the original y,t (theta,x) space.
  ##
  ## The estimated mean E[Y_i | Design ] = mu.h + Ur.h %*% lambda.h^{1/2} %*% EmuMean_i
  ## where i indexes the set of original parameters x1,x2 etc and
  ## lambda.h is the eigenvalue matrix of the pca decomp
  ## ur.h is the eigenvector matrix of same.
  ## mu.h is the mean of the various Y_i curves in the t space.
  ##
  ## The estimated Variance V_ti = V[Y_i | Design] = Sum_L (U_t_l^2 * lambda_l) * V_l_i
  ## where V_ti is the variance of the t'th location (1:ntps) in the pre-pca coordinates and i is indexing
  ## the various pca decomp'ed coordinates (1:nr)
  ##
  ## this can be a little confusing if you end up having as many coord in the emulated x1,x2 space
  ## say 100 as you have locations in the t space.
  ##
  ## 
  ntps <- length(pca.decomp$t)
  nemupts <- dim(emuDesign)[1]
  nr <-  dim(pca.decomp$zmat)[1] # the number of pca cpts we used
  lam <- diag(sqrt(pca.decomp$esys$values[1:nr])) # evalues
  # the reconstructed mean
  yMeanEmuRecon <- pca.decomp$mu  + pca.decomp$ur %*% lam %*% t(zemu.result$emuMean)
  yVarRecon <- matrix(0, nrow=ntps, ncol=nemupts)
  # the reconstructed variance
  for(i in 1:ntps){
    for(j in 1:nemupts){
      #browser()
      # this is a bit confusing we're doing
      # yVar_i_j = Sum_k ( U_i_k * U_i_k * lambda_k * V_k_j)
      yVarRecon[i,j] <- sum((pca.decomp$ur[i,]**2)*(pca.decomp$esys$values[1:nr])*zemu.result$emuVar[j,])
    }
  }
  result <- list(tvalues =pca.decomp$t, mean = yMeanEmuRecon, var = yVarRecon)
  invisible(result)
}



plotReconAll <- function(samples, recon, ylimits=c(0,1)){
  ## pca.decomp: is the pca decomposition of the original design
  ## 
  ## samples: are the curves drawn for various values in the parameter design space, they
  ## have extent along the "t" direction
  ##
  ## reconResult: contains the reconstructed emulated curves, these have extent in the "t" direction
  ## and encode some previously determined variation in the original parameter space.
  plot(recon$tvalues, recon$mean[,1], ylim=ylimits,type="l", xlab="y(rapidity)", ylab="dN_pi/dy")
  grid()
  title(main="Reconstructed Curves and original samples", sub="black are emulated, colored are samples")
    for(i in 2:dim(zemu.points)[1]){
    lines(recon$tvalues, recon$mean[,i])
  }

  for(i in 1:dim(samples$y)[2])
    lines(recon$tvalues, samples$y[,i], col=i, lwd=2)
}


plotReconByIndex <- function(recon, index, first=TRUE, color=1){
  ## just plots the ith result as a mean and variance
  ## if FIRST is true this opens a new plot window, otherwise we just plot in the curent one
  if(first == TRUE){
    plot(recon$t, recon$mean[,index], type="l", xlab="t", ylab="y", lwd=2, col=color)
  } else{
    lines(recon$t, recon$mean[,index], type="l", xlab="t", ylab="y", lwd=2, col=color)
  }
  confUp <- sqrt(1.96*recon$var[,index])+ recon$mean[,index]
  confDown <- recon$mean[,index] - sqrt(1.96*recon$var[,index])
  lines(recon$t, confUp, lty=2, lwd=2, col=2)
  lines(recon$t, confDown, lty=2, lwd=2, col=2)
}

plotAFewRecon <- function(recon){
  ## just draw a few of the curves with their confidence intervals
  plotReconByIndex(recon, 1, TRUE)
  grid()
  plotReconByIndex(recon, 10, color=3, first=FALSE)
  plotReconByIndex(recon, 23, color=4, first=FALSE)
  plotReconByIndex(recon, 50, color=5, first=FALSE)
  plotReconByIndex(recon, 74, color=6, first=FALSE)
}
