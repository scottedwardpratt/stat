#  $Id:$
#
#  Scripts to illustrate PCA emulation approach with toy
#  dataset based on draws from a damped sine wave
#

retval <- require("lhs", quietly=TRUE)                    # Latin Hypercube Sampling pkg
#options(echo=FALSE)
if(retval == FALSE){
  print("installing the lhs library")
  install.packages(c("lhs"))
}

#####################################################
# Test function: damped sine wave
#####################################################
yM <- function(t, x1, x2)
  { return(exp(-outer(t, x1))*sin(outer(t,12*x2))); }

#####################################################
# Function to generate samples
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

#####################################################
# Function to display design
#####################################################
#
# for higher dimensional designs, index1 and 2 let you
# pick which pair to plot against each other
#
view.des <- function(obj=NA, index1=1, index2=2) {
  if(missing(obj)) stop("Compulsary list argument from gen.samp()");
#  x1 <- obj$des[,1];
#  x2 <- obj$des[,2];
  nS <- dim(obj$des)[2];
  plot(t(obj$des[index1,]), t(obj$des[index2,]), xlab=expression(x[index1]), ylab=expression(x[index2]));
  title(bquote(paste( .(nS), " Design points")));
}
  
#####################################################
# Function to find PCA
#####################################################
emu.pca <- function(sample=NA, nr=3, ifile="", silent=TRUE) {
  if(missing(sample)) {
    if(missing(ifile))
      stop("Must give sample OR ifile")
    else
      load(ifile);
  }
  # One way or another we should have list "sample" now
  des   <- sample$des;        # Dim: 2  x nS
  ts    <- sample$t;          # Dim: nS x  1 (no! it's ntps)
  ys    <- sample$y;          # Dim: ntps x 
  ntps  <- dim(ys)[1];
  nreps <- dim(ys)[2];
  mu.h  <- apply(ys, 1, mean);
  sig.h <- (ys-mu.h) %*% t(ys-mu.h) / nreps;
  esys  <- eigen(sig.h);
  ur.h  <- esys$vec[,1:nr];
  lam.h <- esys$val[1:nr];
  lam.inv.sqrt <- diag(1/sqrt(lam.h));

  if(!silent){
  # Demonstrate how choice of rv holds up:
    top.r <- min(nr + 5, ntps);
    frac  <- cumsum(esys$val[1:top.r]) / sum(esys$val);
    prv   <- plot(1:top.r, frac, xlab="Index", ylab="Eigenvalue",type="b");
    abline(v=nr,col="blue");
    abline(h= 1,col="red")
    title(bquote(paste("Top ",.(nr)," evs are ", .(round(100*frac[nr],4)),
                       "% of total"))); #,cex=1.2);
    cat("Press Enter to continue...");
    readline();
    plot(ts, ur.h[,1], col=1, type="l", xlab="t", ylab="U", ylim=range(ur.h));
    for(i in 2:nr)
      lines(ts, ur.h[,i], col=i);
    title(main = "showing the first few principle cpts")
  }

  zmat  <- lam.inv.sqrt  %*% t(ur.h) %*% (ys - mu.h);
  rv    <- list(des=des, t=ts, y=ys, nreps=nreps, ntps=ntps, mu=mu.h,
                sig=sig.h, esys=esys, ur.h=ur.h, zmat=zmat
             );
  invisible(rv);
}

## creates a lovely grid
# thank you expand grid, this is the best solution i've seen to this problem
createEmuPoints <- function(ndim, npointsSide, rangemin=0, rangemax=1){
  if(ndim == 1){
    ans <- seq(from=rangemin, to=rangemax,length.out=npoints)
  } else if(ndim >1){
    side <- list(c(0:(npointsSide-1)))
    ans <- ((rangemax-rangemin)/npointsSide)*as.matrix(expand.grid(rep(side, ndim)))+rangemin
  }
}


estimateThetas <- function(pca.decomp){
# estimates the hyperparameters for the emulator
# the return is a matrix whose rows are the hyperparameters for each cpt of the z decomposition (there should be nr rows therefore)
# the args are:
#
# design: a matrix of the points at which the training values (zmat) are evaluated, these are the various params that the
# model was run at. A matrix npts * nparams
# npts: the number of points in the parameter space of the model.
# training values: a matrix of the values of the model (the zcpts) at each element of the design. npts * nr
# the number of training values: nr
  
  thetas.est <- multidim(t(pca.decomp$des), pca.decomp$nreps, t(pca.decomp$zmat), dim(pca.decomp$zmat)[1])
  # we'll save the thetas because doing the estimation rapidly gets boring
  #write.table(thetas.est, "thetas.dat")
  invisible(thetas.est)
}

reconCurveAtPoint <- function(point, thetas, pca.decomp){
  nr <- dim(pca.decomp$zmat)[1]
  nparams <- dim(pca.decomp$des)[1]
  nmodelpts <- dim(pca.decomp$des)[2]
  emuResult <- list(xpos = point, mean = rep(0, nr), var = rep(0, nr))
  ntps <- length(pca.decomp$t)
  for(i in 1:nr){
    model <- list(xmodel=as.matrix(t(pca.decomp$des)), training=t(pca.decomp$zmat)[,1])
    temp <- callEmulateAtPoint(model, as.numeric(thetas[i,]), as.numeric(point), nmodelpts, nparams=nparams, nthetas=(nparams+2))
    emuResult$mean[i] <- temp$mean
    emuResult$var[i] <- temp$var
  }
  # now we need to recon the curve
  yMeanEmuRecon <- pca.decomp$mu + pca.decomp$ur %*% diag(sqrt(pca.decomp$esys$values[1:nr])) %*%
    emuResult$mean

  yVarEmuRecon <- rep(0, ntps)
  
  for(i in 1:ntps){
    ## ccs
    ## not sure this is right, rederive
    yVarEmuRecon[i] <- sum((pca.decomp$ur[i,]**2)*(pca.decomp$esys$values[1:nr])*emuResult$var)
  }

  finalResult <- list(xpos = point, tvalues = pca.decomp$t, mean=yMeanEmuRecon, var = yVarEmuRecon)  
  invisible(finalResult)
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


