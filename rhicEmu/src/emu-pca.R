#  $Id:$
#
#  Scripts to illustrate PCA emulation approach with toy
#  dataset based on draws from a damped sine wave
#

retval <- require("lhs", quietly=TRUE)                    # Latin Hypercube Sampling pkg
#options(echo=FALSE)
if(retval == FALSE){
  print("#installing the lhs library")
  install.packages(c("lhs"))
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
#
# ccs, we should add or modify this to dynamically
# decide the minimum number of eigenvalues to get above a certain
# threshhold.
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


##
## project the emulators evaluated at point back into the functional
## space using the pca.decomp calculated above
##
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
    ## this is suspect
    ## 
    yVarEmuRecon[i] <- sum((pca.decomp$ur[i,]**2)*(pca.decomp$esys$values[1:nr])*emuResult$var)
  }

  finalResult <- list(xpos = point, tvalues = pca.decomp$t, mean=yMeanEmuRecon, var = yVarEmuRecon)  
  invisible(finalResult)
}



## same as reconCurveAtPoint but for some set of points
##
## the results back from the C code are stored in the list
## emuResults, each column is the result from one of the z-emulators at each of
## the nEmupts
##
## need to compare this agains reconCurveAtPoint, do they agree?
reconCurveAtList <- function(pointSet, thetas, pca.decomp){
  nr <- dim(pca.decomp$zmat)[1]
  nparams <- dim(pca.decomp$des)[1]
  nmodelpts <- dim(pca.decomp$des)[2]
  # points = nemupts x nparams points
  nemupts <- dim(pointSet)[1]
  emuResult <- list(xpos = pointSet, mean = matrix(0, ncol=nr, nrow=nemupts), var=matrix(0, ncol=nr, nrow=nemupts))
  ntps <- length(pca.decomp$t)
  
  #print(emuResult$mean)
  
  for(i in 1:nr){
    model <- list(xmodel=as.matrix(t(pca.decomp$des)), training=t(pca.decomp$zmat)[,1])
    
    temp <- callEmulateAtList(model, as.numeric(thetas[i,]), as.numeric(pointSet), nemupts, nmodelpts, nparams=nparams, nthetas=(nparams+2))
    #print(temp)
    # these now come back as vectors
    ## print(length(temp$mean))
    ## print(length(temp$var))
    ## print(dim(emuResult$mean))
    ## print(dim(emuResult$var))
    emuResult$mean[,i]<- temp$mean   
    emuResult$var[,i] <- temp$var
  }
  #print(emuResult$mean)
  # now we need to recon the curve
  #
  # we read each result as a row of the emuResult$Mean matrix and
  # multiply it into the yMeanEmuRecon matrix, which is stored with each curve
  # in a single column
  #browser()
  yMeanEmuRecon <- matrix(0, ncol=nemupts, nrow=ntps)
  for(i in 1:nemupts){
    yMeanEmuRecon[,i] <- pca.decomp$mu + pca.decomp$ur %*% diag(sqrt(pca.decomp$esys$values[1:nr])) %*%
      emuResult$mean[i,]
  }

  #print(yMeanEmuRecon)
  
  yVarEmuRecon <- matrix(0, ncol=nemupts, nrow=ntps)

  for(i in 1:ntps){
    for(j in 1:nemupts){
      #browser()
      # this is a bit confusing we're doing
      # yVar_i_j = Sum_k ( U_i_k * U_i_k * lambda_k * V_k_j)
      yVarEmuRecon[i,j] <- sum((pca.decomp$ur[i,]**2)*(pca.decomp$esys$values[1:nr])*emuResult$var[j,])
    }
  }

  finalResult <- list(xpos = pointSet, tvalues = pca.decomp$t, mean=yMeanEmuRecon, var = yVarEmuRecon)  
  invisible(finalResult)

}


##
## generate the emulator data without projecting back
##
## similar to recon curve at list, but does not actually
## do the projection back into the functional space
##
emulateAtListNoProject <- function(pointSet, thetas, pca.decomp){
  nr <- dim(pca.decomp$zmat)[1]
  nparams <- dim(pca.decomp$des)[1]
  nmodelpts <- dim(pca.decomp$des)[2]
  # points = nemupts x nparams points
  nemupts <- dim(pointSet)[1]
  emuResult <- list(xpos = pointSet, mean = matrix(0, ncol=nr, nrow=nemupts), var=matrix(0, ncol=nr, nrow=nemupts))
  ntps <- length(pca.decomp$t)
  
  #print(emuResult$mean)
  
  for(i in 1:nr){
    model <- list(xmodel=as.matrix(t(pca.decomp$des)), training=t(pca.decomp$zmat)[,1])
    
    temp <- callEmulateAtList(model, as.numeric(thetas[i,]), as.numeric(pointSet), nemupts, nmodelpts, nparams=nparams, nthetas=(nparams+2))
    #print(temp)
    # these now come back as vectors
    ## print(length(temp$mean))
    ## print(length(temp$var))
    ## print(dim(emuResult$mean))
    ## print(dim(emuResult$var))
    emuResult$mean[,i]<- temp$mean   
    emuResult$var[,i] <- temp$var
  }

  finalResult <- list(xpos = pointSet, mean=emuResult$mean, var = emuResult$var)  
  invisible(finalResult)
}

