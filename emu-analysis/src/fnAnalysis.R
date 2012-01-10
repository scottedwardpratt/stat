##
## C.C-S, cec24@phy.duke.edu
## sep 2011
##
## functions to construct and sample an emulator based upon some
## P observables which are not a-priori independent
##
## we assuming no particular ordering in the functional (t) space
## (not necessarily some kind of spectrum etc)
## 
##
## 0) arrange data into a suitable rep
## 
## 1) construct a pca decomp of the data, reducing from P to r
##    independent components.
## 
## 2) construct an emulator for each of the components
##
## 3) compute the mean / var for each of the original observables,
## validate the constructed emulator against with-held data
## 
## 4) compute the combined implausibility across the design space,
## project this in various ways...
##
##
## step by step to get to this level of started we need to:
## 1) define at the global scope expData, modelData, designData objects
##
## expData -> list((vector)obsValue, (vector)obsError), this defines the distributions of the observables we are going
## to develop a feasible region for. obsValue is the reported mean, obsError should be the s.d. Currently
## we lack support for other distributions, but we can deal with this when we generate the fn density.observable
##
## modelData -> matrix(nrows=nmodelpts, ncols=nbins), contains the sampled values of our model over the design,
## each row represents a single sample of the nbins-variate output.
##
## designData -> matrix(nrows=nmodelpts, ncols=nparams), the locations in the design space where
## we generated modelData
## 
## we also need the following: nparams (dimensionality of the design space), nbins (number of observables we're treating),
## desNames (a vec of the names of the design variables, for plots), 
##
## 2) generate a sample, the return will be the basis of our functional data object that we'll keep
## passing around (its our state) so:
## fnData <- fn.sample.gen(cov.fn=1, reg.order=1) # choice of cov.fn and reg.order get kept forever...
##
## 3) pca.decompose the sample, you can set the variance cutoff (high values are suggested for small nbins) this
## will determine how many eigenvalues (dimensions) are thrown away
## fnData.pca <- fn.pca.gen(fnData, cutOff=0.98)
##
## 4) estimate the hyperparams for this fnData model 
## fnData.final <- fn.estimate(fnData.pca)
##
## 5) now emulate the model at some locations and then compute the implausibility  etc
## eg:
## stepData <- fn.emu.steps(fnData.final, ...)
## impSteps <- fn.implaus.steps(fnData.final, stepData)
## 


## arranging the data
## make this compatible with emuOver and impOver so we'll
## suppose that the following items exist at the global level
## 
## - nbins: the number of observables 
## - nruns:  how many samples we have
## - modelData matrix(nrows=nruns, ncols=nbins)
## - expData list(obsValue(vector length nbins), obsError(vector length nbins))
## - designData
##
## @param cov.fn.in we'll use the same cov fn throughout, so it should only be
## set here
## @param reg.order.in set the regression order for the whole analysis
##
## @return list(model.sample, exp.obs, cov.fn, reg.order, nruns, nbins)
fn.sample.gen <- function(cov.fn.in=1, reg.order.in=1 ){
  ## we'll scale the model.data and design, currently
  ## do this using the R builtin
  ## modelData.scaled <- scale.unit(modelData)
  ## des.scaled <- scale.unit(designData)
  modelData.scaled <- scale(modelData)
  des.scaled <- scale(designData)

  ## the classic emu.pca object
  sample <- list(t=seq(1,nbins), des=t(des.scaled), y=t(modelData.scaled))

  retVal <- list(model.sample=sample, exp.obs=expData, cov.fn=cov.fn.in, reg.order=reg.order.in, nruns=nruns, nbins=nbins, nparams=dim(designData)[2])

  invisible(retVal)
}





## takes data.in and scales its columns uniformly so that
## on return the values fall within a unit hypercube
##
## @param data.in is a matrix, and we'll scale it column by column
scale.unit <- function(data.in){
  nbins <- dim(data.in)[2]
  nrows <- dim(data.in)[1]
  scales <- rep(NA, nbins)
  centers <- rep(NA, nbins)
  data.scaled <- matrix(0, nrow=nrows, ncol=nbins)
  for(i in 1:nbins){
    scales[i] <- max(data.in[,i]) - min(data.in[,i])
    centers[i] <- min(data.in[,i])
    data.scaled[,i] <- (data.in[,i] - centers[i])/scales[i]
  }
  attr(data.scaled, "scaled:scale") <- scales
  attr(data.scaled, "scaled:center") <- centers

  invisible(data.scaled)
}
  


##
## create the pca.decomposition using emu.pca.adapt
##
## @param fn.data, a list created by fn.sample.gen
## @return fn.data, same as arg but with the pca.decomposition object added
## 
##
## pca.decomp    <- list(des=des, t=ts, y=ys, nreps=nreps, ntps=ntps, mu=mu.h,
##              sig=sig.h, esys=esys, ur.h=ur.h, zmat=zmat )
##
## pca.decomp$des experimental design (id to fn.data$des)
## pca.decomp$ts functional space index set (id to fn.data$model.sample$t)
## pca.decomp$ys fn.data$model.sample$y
## pca.decomp$nreps number of realisations of the fnal data (nruns)
## pca.decomp$ntps number of sites in the functional space (nbins)
## pca.decomp$mu.h sample mean of Ys
## pca.decomp$sig.h sample variance of Ys
## pca.decomp$esys eigensystem of pca decomposition (important)
## pca.decomp$ur eigenvectors of pca.decomp
## pca.decomp$zmat matrix of training values for the r components
## 
## 
fn.pca.gen <- function(fn.data, cutOff=0.96, nrGuess=3,silent.in=TRUE ){
  ## 
  pca.result <- emu.pca.adapt(sample=fn.data$model.sample, cutOff=cutOff,
                              nrGuess=nrGuess, silent=silent.in)


  fn.data$pca.decomp <- pca.result
  ## return the list, now containing the pca.decomposition
  invisible(fn.data)
}

##
## estimate the hyperparams for our model and store them in
## fn.data
##
## @return fn.data contains the element thetas.est, a matrix of the
## estimated hyperparameters for the emulator, each row is a separate
## component, each cpt must use the same cov.fn & reg.model (for now)
##
fn.estimate <- function(fn.data){
  if(is.null(fn.data$pca.decomp)){
    stop("fn.estimate needs a pca decomposition to estimate")
  }

  # the transposes are important (annoyingly)
  thetas <- multidim(model=t(fn.data$model.sample$des),
                     training=t(fn.data$pca.decomp$zmat),
                     nmodelpts=fn.data$nruns, nydims=fn.data$pca.decomp$nr,
                     cov.fn=fn.data$cov.fn, reg.order=fn.data$reg.order)

  fn.data$thetas.est <- thetas

  invisible(fn.data)
}


## validate the given emulator (fn.data) by comparing predictions at holdback.data$des
## with their true values holdback.data$training
##
## we'll use the Mahalonobis distance as in implaus.joint
fn.validate <- function(fn.data, holdback.data){

  nvalpts <- dim(holdback.data$des)[1]
  nparams <- dim(holdback.data$des)[2]
  nobs <- dim(holdback.data$training)[2]
  
  ## for unscaling the mean
  trainCenter.vec <- attr(fn.data$model.sample$y, "scaled:center")
  trainScale.vec <- attr(fn.data$model.sample$y, "scaled:scale")

  ## for unscaling the variance
  vScale.mat <- matrix(0, nrow=nbins, ncol=nbins)
  for(i in 1:nbins)
    for(j in 1:nbins)
      vScale.mat[i,j] <- trainScale.vec[i] * trainScale.vec[j]

  ## to scale the test-design points
  desCenter.vec <- attr(fn.data$model.sample$des, "scaled:center")
  desScale.vec <- attr(fn.data$model.sample$des, "scaled:scale")

  test.des <- holdback.data$des
  ## undo the normal scaling on the design first
  for(i in 1:nparams){
    test.des[,i] <- (test.des[,i] - desCenter.vec[i]) / desScale.vec[i]
  }

  ## generate predictions
  ## 
  ## clist = list(xpos, tvalues, mean, var)
  ## where xpos == pointList, tvalues == pca.decomp$t
  clist <- reconCurveAtList(test.des, fn.data$thetas.est, fn.data$pca.decomp,
                            cov.fn=fn.data$cov.fn, reg.order=fn.data$reg.order)

  deviance.joint <- rep(NA, nvalpts)
  deviance.single <- matrix(0, nrow=nvalpts, ncol=nobs)
  
  for(i in 1:nvalpts){
    V.mat <- (vScale.mat)*clist$tCovList[[i]]
    V.mat.inv <- solve(V.mat)
    mean.vec <- clist$mean[,i] * trainScale.vec + trainCenter.vec
    deviance.joint[i] <- t(holdback.data$training[i,] - mean.vec) %*% V.mat.inv %*% (holdback.data$training[i,] - mean.vec)
    for(j in 1:nobs){
      deviance.single[i, j] <- (holdback.data$training[i,j] - mean.vec[j]) ** 2 / (V.mat[j,j])
    }
  }
  
  result <-  list(joint=deviance.joint, single=deviance.single)
}

## compute the emulated mean and var for the reconstructed
## functional variables over a grid in dimA, dimB in design space
## a list of covaraince matrices between the different functional
## observables is also saved in tCovList
##
## this is all done in the scaled variables
## 
## @param fn.data, must contain thetas.est and pca.decomp
## @param dimA,dimB dimenstions in the design space to emulate over
## @param fixedVals either NULL (if nparams==2) or a vector of nparams length, containing
## NA's at the dimA and dimB slots.
## 
## @param nemupts number of grid sites per side
##
## @return list(point_list, mean, var, tCovList, cov.fn, reg.order, obs.dim, dimA, dimB)
## 
fn.emulate.slice <- function(fn.data, dimA, dimB, fixedValVec=NULL, nemupts=32){
  if(is.null(fn.data$pca.decomp) || is.null(fn.data$thetas.est)){
    stop("fn.emulate needs a pca.decomp && a set of estimated thetas")
  }

  # recall that the model is transposed in the
  # fn.data
  des.scaled <- t(fn.data$model.sample$des)
  ## now create the correct ranges
  rangeA <- c(min(des.scaled[,dimA]), max(des.scaled[,dimA]))
  rangeB <- c(min(des.scaled[,dimB]), max(des.scaled[,dimB]))
  stepSizeA <- (rangeA[2] - rangeA[1]) / (nemupts)
  stepSizeB <- (rangeA[2] - rangeA[1]) / (nemupts)

  pointList <- matrix(0, nrow=nemupts**2, ncol=nparams)

  for(i in 1:nemupts){
    for(j in 1:nemupts){
      pointVec <- rep(NA, nparams)
      pointVec[dimA] <- rangeA[1]+i*stepSizeA
      pointVec[dimB] <- rangeB[1]+j*stepSizeB

      ## copy in the rest of the non na entries from
      ## fixedValVec 
      for(vecIndex in 1:nparams){
        if( is.na(pointVec[vecIndex]) == TRUE){
          pointVec[vecIndex] <- fixedValVec[vecIndex]
        }
      }

      pointList[j+nemupts*(i-1),] <- pointVec
    }
  }
  
  ## clist = list(xpos, tvalues, mean, var)
  ## where xpos == pointList, tvalues == pca.decomp$t
  clist <- reconCurveAtList(pointList, fn.data$thetas.est, fn.data$pca.decomp,
                            cov.fn=fn.data$cov.fn, reg.order=fn.data$reg.order)

  ## tCovList is the list of t-space covariance matrices,
  ## these are a posterior estimate on the covariance between the
  ## different observables in the functional space
  result <- list(pointList=pointList,
                 mean=clist$mean, 
                 var=clist$var,
                 tCovList=clist$tCovList,
                 reg.order=fn.data$reg.order, cov.fn=fn.data$cov.fn,
                 dimA=dimA, dimB=dimB, fixedVals=fixedValVec,
                 npts=nemupts)

  invisible(result)
}


## sample the density of our multivariate emulator at a given location in the parameter space
## we'll need this for estimating the functional region
##
## currently the results from the emulator are not rescaled
##
## we need the library "mvtnorm" installed
##
## this runs using the multivariate extension to the quick emulate.MC
## functions developed for the post-sample pilot project
##
## @param sample.pt the location in the design space we want to evaluate the density at
## @param sample.value the vector (nydims) long of observables we want to compute
## the density of.
##
## @requires: an already estimated fn.data sample
##
## need to be careful about transposing the data
fn.sampleEmu.density <- function(fn.data, sample.pt, sample.value, use.log=TRUE){

  pca.decomp <- fn.data$pca.decomp
  

  ## we need to rotate the mean and compute the covariance matrix
  ##
  ## shortcut some numbers we'll use
  nr <- dim(pca.decomp$zmat)[1]
  nparams <- dim(pca.decomp$des)[1]
  nmodelpts <- dim(pca.decomp$des)[2]
  ntps <- length(pca.decomp$t)

  ## emulate a single point, this result in the functional space!
  ## this needs to return a mean and var vector
  result <- emulate.MC.multi(sample.pt, nr)

  yMeanEmuRecon <- rep(0, ntps) ## we have ntps different observables to reconstruct into
  yVarMatEmuRecon <- matrix(0, nrow=ntps, ncol=ntps) ## the full covar matrix at sample.pt

  yMeanEmuRecon <- pca.decomp$mu + pca.decomp$ur %*% diag(sqrt(pca.decomp$esys$values[1:nr])) %*% result$mean

  for(i in 1:ntps){
    for(j in 1:ntps){
      yVarMatEmuRecon[i,j] <- sum(pca.decomp$ur[i,]*pca.decomp$ur[j,]*
                                  pca.decomp$esys$values[1:nr] * result$var)
    }
  }

  density <- dmvnorm(sample.value, mean=yMeanEmuRecon,
                     sigma=yVarMatEmuRecon, log=use.log)
  invisible(density)
}



## compute the joint implausibility for a slice
## and also the independent implaus for each obs
##
## implaus_joint^2 = (E(f(x) - z))^t V^{-1} E(f(x) - z)
## where E(f(x)) is a vector of the emulator mean at the point x
## in the design space.
## 
## V is the matrix of space posterior covariances plus any variances associated
## with the model and real-observations along the diagonal...
## 
## uses the globally available expData (hopefully)
fn.implaus.slice <- function(fn.data, slice.data){
  ## for unscaling the mean, in slice.data
  trainCenter.vec <- attr(fn.data$model.sample$y, "scaled:center")
  trainScale.vec <- attr(fn.data$model.sample$y, "scaled:scale")

  ## for unscaling the variance
  vScale.mat <- matrix(0, nrow=nbins, ncol=nbins)
  for(i in 1:nbins)
    for(j in 1:nbins)
      vScale.mat[i,j] <- trainScale.vec[i] * trainScale.vec[j]
      
  
  npts <- slice.data$npts

  implaus.joint <- matrix(0, nrow=npts, ncol=npts)
  implaus.inde.max <- matrix(0, nrow=npts, ncol=npts)
  implaus.inde <- vector("list", fnData$nbins)
  for(i in 1:fnData$nbins)
    implaus.inde[[i]] <- matrix(0, nrow=npts, ncol=npts)

  for(i in 1:npts){
    for(j in 1:npts){
      listIndex <- j+npts*(i-1)
      V.mat <- (vScale.mat)*slice.data$tCovList[[listIndex]] + diag(expData$obsError**2)
      V.mat.inv <- solve(V.mat)
      ## remove the normal scaling before undoing the centering
      mean.vec <- slice.data$mean[,listIndex] * trainScale.vec + trainCenter.vec
      implaus.joint[i,j] <- t(expData$obsValue - mean.vec) %*% V.mat.inv %*% (expData$obsValue - mean.vec) 
      imp.inde.vec <- NULL
      for(k in 1:fnData$nbins){
        implaus.inde[[k]][i,j] <- (mean.vec[k] - expData$obsValue[k])**2 / (V.mat[k,k]) 
        if(implaus.inde[[k]][i,j] > 4.0)
          implaus.inde[[k]][i,j] <- 4.0
        imp.inde.vec <- c(imp.inde.vec, implaus.inde[[k]][i,j])
      }
      implaus.inde.max[i,j] <- max(imp.inde.vec)
    }
  }

  result <- list(implaus.joint=implaus.joint,
                 implaus.inde=implaus.inde,
                 implaus.inde.max=t(implaus.inde.max),
 npts=npts)
  
}

## compute the implausibility over a set of slices through the model
## the slices (stepData) have already been created by a call to fn.emu.steps
##
## @param stepData the set of slices we want to compute the implaus on
## @return a list with the implaus on for each slice
fn.implaus.steps <- function(fn.data, stepData){
  nslices <- stepData$nsteps

  imp.list <- vector("list", nslices)
  for( i in 1:nslices){
    imp.list[[i]] <- fn.implaus.slice(fn.data, stepData[[i]])
  }

  imp.list$npts <- imp.list[[1]]$npts
  imp.list$nsteps <- stepData$nsteps
  imp.list$stepDim <- stepData$stepDim
  imp.list$dimA <- stepData$dimA
  imp.list$dimB <- stepData$dimB
  imp.list$fixedVals <- stepData$fixedValVec
  imp.list$stepVec <- stepData$stepVec

  
  invisible(imp.list)
}



## plot an emulated slice produced by fn.emulate.slice for obsDim
##
## by default we unscale the data back to its real size
##
## @param plot.thing is a vector of nemupts**2 which is plotted
## over designA, designB
fn.plot.slice <- function(fn.data, slice.data, plot.thing, unscale=TRUE, legend.in="",
                          xlabel="", ylabel="", plot.des=FALSE, imp=FALSE){

  designA <- fn.data$model.sample$des[slice.data$dimA,]
  designB <- fn.data$model.sample$des[slice.data$dimB,]
  
  if(imp==FALSE){
    plot.data <- matrix(plot.thing, nrow=slice.data$npts, ncol=slice.data$npts, byrow=TRUE)
  } else {
    plot.data = plot.thing
  }
  #mean <- matrix(slice.data$mean[obsIndex,], nrow=slice.data$npts, ncol=slice.data$npts)
  #var <- matrix(slice.data$var[obsIndex,], nrow=slice.data$npts, ncol=slice.data$npts)

  range1 <- seq(min(designA), max(designA), length=slice.data$npts)
  range2 <- seq(min(designB), max(designB), length=slice.data$npts)

  if(unscale==TRUE){ # we'll return to the original scale for plotting

    desACenter <- attr(fn.data$model.sample$des, "scaled:center")[slice.data$dimA]
    desBCenter <- attr(fn.data$model.sample$des, "scaled:center")[slice.data$dimB]

    desAScale <- attr(fn.data$model.sample$des, "scaled:scale")[slice.data$dimA]
    desBScale <- attr(fn.data$model.sample$des, "scaled:scale")[slice.data$dimB]

    range1 <- range1 * desAScale + desACenter
    range2 <- range2 * desBScale + desBCenter

    designA <- designA * desAScale + desACenter
    designB <- designB * desBScale + desBCenter

  }

  if(imp==FALSE){
    image(range1, range2, plot.data, axes=FALSE, col=heat.colors(16), xlab=xlabel, ylab=ylabel)
    contour(range1, range2, plot.data, nlevels=10, col="black", add=TRUE, cex.lab=0.5, labcex=0.8)
  } else {
    if(max(plot.data) <= 4.0){
      ## plot the implausibility on a nice set of breaks
      ils <- c(0.1, 0.2, 0.5, 1.0, 1.5, 2.0)
    } else {
      ils <- c(0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0, 16.0)
    }
    breaks <- c(0, ils)
    nlevels <- length(ils)
    image(range1, range2, plot.data, axes=FALSE, col=rev(heat.colors(nlevels)), breaks=breaks,
          xlab=xlabel, ylab=ylabel)
    contour(range1, range2, plot.data, levels=ils, labels=ils, col="black", add=TRUE, cex.lab=0.5, labcex=0.8)
  }

  legend("topright", legend.in, bg="white", cex=1.5, bty="n")
  if(plot.des==TRUE){
    points(designA, designB, pch=3)
    title(xlab=xlabel, ylab=ylabel, outer=TRUE, cex.lab=3.0)
  }
  axis(1, cex.axis=1.5)
  axis(2, cex.axis=1.5)

}

## 
## emulate a set of in planes of dimA and dimB, stepping through stepDim
## the rest of the parameters are held at their fixedValVec values.
##
## @param fixedValVec should be a vector of length nparams, entries which
## are to be varied over (dimA, dimB, stepDim) should be NA, the remaining
## entries should be held at their fixed values
##
## @param dimA, dimB index of the dimensions to emulate the model over
## @param stepDim, index of the dimension to march through the space
## @param hyperSLice, if true each point in the slice is generated through sampling in the directions orthogonal
## to the cube we're building
fn.emu.steps <- function(fn.data, dimA, dimB, stepDim, fixedValVec=NULL, nsteps=9,
                         range.Min=NULL, range.Max=NULL, hyperSlice=FALSE, nemupts=32){

  if(is.null(range.Min) && is.null(range.Max)){
    minVal <- min(fn.data$model.sample$des[stepDim,])
    maxVal <- max(fn.data$model.sample$des[stepDim,])
  } else {
    minVal <- range.Min
    maxVal <- range.Max
  }

  stepSize <- (maxVal - minVal) / nsteps

  sliceList <- vector("list", nsteps)

  # store the values of the stepping parameter
  stepVec <- rep(NA, nsteps)
  
  for(i in 1:nsteps){
    stepVal <- minVal + stepSize * (i-1)
    ## pass this to fn.emulate.slice, remaining NA fields should be at dimA and dimB
    fixedStep <- rep(NA, nparams)
    ## there are some (nparams  - 3) entries in fixedValVec
    ## we replace the stepDim entry with the value its fixed at for this
    ## sweep in 2d
    fixedStep[stepDim] <- stepVal
    if(is.null(fixedValVec) == FALSE){
      # copy in the rest of the fixed values, the non.na values in 
      # fixedValVec
      for(fixedIndex in 1:nparams){
        if(is.na(fixedValVec[fixedIndex]) == FALSE){
          fixedStep[fixedIndex] <- fixedValVec[fixedIndex]
        }
      }
    }
    if(hyperSlice == FALSE){
      sliceList[[i]] <- fn.emulate.slice(fn.data, dimA, dimB, fixedValVec=fixedStep, nemupts=nemupts)
    } else {
      sliceList[[i]] <- fn.emulate.hyper.slice(fn.data, dimA, dimB, stepDim, stepVal)
    }
    stepVec[i] <- stepVal
  }
  
  sliceList$nsteps <- nsteps
  sliceList$stepDim <- stepDim
  sliceList$dimA <- dimA
  sliceList$dimB <- dimB
  sliceList$fixedVals <- fixedValVec
  sliceList$stepVec <- stepVec
  # now return the list of slices we created
  invisible(sliceList)
}

## plot the mean or var of a set of slices
## through the functional data for a given
## observable
fn.plot.steps <- function(fn.data, step.list, obsIndex, plot.var=FALSE){
  par(mfrow=c(3,3), mar=c(1,1,0,0), oma=c(4,5,3,0))


  trainCenter <- attr(fn.data$model.sample$y, "scaled:center")[obsIndex]
  trainScale <- attr(fn.data$model.sample$y, "scaled:scale")[obsIndex]

  for(i in 1:step.list$nsteps){
    ## undo the normal scaling here
    stepValue <- step.list$stepVec[i]
    stepValue.unscaled <- stepValue * attr(fn.data$model.sample$des, "scaled:scale")[step.list$stepDim] +
      attr(fn.data$model.sample$des, "scaled:center")[step.list$stepDim]
      
    leg.buffer <- paste(desNames[step.list$stepDim], ": ", round(stepValue.unscaled,2), sep="")

    mean.scaled <- step.list[[i]]$mean[obsIndex,] * trainScale + trainCenter
    
    if(plot.var==TRUE){
      var.scaled <- step.list[[i]]$var[obsIndex,] * (trainScale**2)
      plot.data <- var.scaled
    } else {
      plot.data <- mean.scaled
    }
    
    if(i==7){
      
      fn.plot.slice(fn.data, step.list[[i]], plot.thing=plot.data, 
                    xlabel=desNames[step.list$dimA], ylabel=desNames[step.list$dimB],
                    legend.in=leg.buffer, plot.des=TRUE)
    } else {
      fn.plot.slice(fn.data, step.list[[i]], plot.thing=plot.data, 
                    xlabel="", ylabel="",
                    legend.in=leg.buffer)
    }
  }
  par(mfrow=c(1,1))
  if(is.null(step.list$fixedVals)==FALSE){
    if(plot.var==FALSE){
      buffer <- paste("E(obs: ", obsIndex, ") fixed: ", step.list$fixedVals)
    } else {
      buffer <- paste("Var(obs: ", obsIndex, ")fixed: ", step.list$fixedVals)
    }
  } else {
    buffer <- paste("E(obs: ", obsIndex, ")")
  }
  title(buffer, outer=TRUE)

}


fn.plot.imp.steps <- function(fn.data, step.list, obsIndex, plot.joint=TRUE, title.in=NULL){
  par(mfrow=c(3,3), mar=c(3,3,0,0), oma=c(4,4,5,0))


  for(i in 1:step.list$nsteps){
    ##
    stepValue <- step.list$stepVec[i]
    stepValue.unscaled <- stepValue * attr(fn.data$model.sample$des, "scaled:scale")[step.list$stepDim] +
      attr(fn.data$model.sample$des, "scaled:center")[step.list$stepDim]
      
    leg.buffer <- paste(desNames[step.list$stepDim], ": ", round(stepValue.unscaled,2), sep="")
    
    if(plot.joint==TRUE){
      plot.data <- step.list[[i]]$implaus.joint
    } else {
      plot.data <- step.list[[i]]$implaus.inde[[obsIndex]]
    }
    
    if(i==7){
      fn.plot.slice(fn.data, step.list, plot.thing=plot.data, 
                    xlabel=desNames[step.list$dimA], ylabel=desNames[step.list$dimB],
                    legend.in=leg.buffer, plot.des=TRUE, imp=TRUE)
    } else {
      fn.plot.slice(fn.data, step.list, plot.thing=plot.data, 
                    xlabel="", ylabel="",
                    legend.in=leg.buffer, imp=TRUE)
    }
  }
  par(mfrow=c(1,1))

  if(plot.joint == FALSE){
    if(is.null(step.list$fixedVals)==FALSE){
      buffer <- paste(title.in, "obs: ", obsIndex, "fixed: ", fixedVals)
    } else {
      buffer <- paste(title.in, "obs: ", obsIndex)
    }
  } else {
    buffer <- paste(title.in, "joint implaus")
  }
 
  title(buffer, outer=TRUE, cex=1.1)

}




#####################################################
##
## functions after this should not be directly called
##
## is there good module support in R?
## 
#####################################################

#####################################################
# Function to find PCA with a variable number of
# components
#
# ccs, madai-meeting-2
#####################################################
# compute a pca decomposition with an adaptive number
# of principle components
#
# supply either a sample object or a file to find one,
# cutOff: the relative fraction of "eigen-power" we
# want our Nr components to find
#
# nrGuess: the minimum number of pca cpts (we think)
# should be important
emu.pca.adapt <- function(sample=NA, cutOff=0.95, nrGuess=3,
                          ifile="", silent=TRUE)
{
  if(missing(sample)){
    if(missing(ifile))
      stop("emu.pca.adapt: give either sample or ifile")
    else
      load(ifile)
  }

  des <- sample$des
  ts <- sample$t
  ys <- sample$y
  ntps <- dim(ys)[1]
  nreps <- dim(ys)[2]
  # the sample mean
  mu.h <- apply(ys, 1, mean)
  # the sample variance
  sig.h <- (ys - mu.h) %*% t(ys - mu.h) / nreps

  # eigensystem
  esys <- eigen(sig.h)
  evecs <- esys$vec
  evals <- esys$val

  frac <- 0
  totalVar <- sum(esys$val)
  count <- 0
  overFlow <- nreps
    
  while(frac < cutOff || count > overFlow){
    frac <- sum(esys$val[1:nrGuess]) / totalVar
    nrGuess <- nrGuess + 1
    if(!silent){
      cat("# Frac: ", frac , " nrGuess: ", nrGuess, "\n")
    }
  }
  nrFinal <- nrGuess
  
  cat("# Final Nr: ", nrFinal, "\n")
  cat("# Final Variance fraction: ", frac, "\n")

  # now construt what we need to match the results of emu.pca
  ur.h <- esys$vec[,1:nrFinal]
  lam.h <- esys$val[1:nrFinal]
  lam.inv.sqrt <- diag(1/sqrt(lam.h))

  if(!silent){
  # Demonstrate how choice of rv holds up:
    top.r <- min(nrFinal + 5, ntps);
    frac  <- cumsum(esys$val[1:top.r]) / sum(esys$val);
    prv   <- plot(1:top.r, frac, xlab="Index", ylab="Eigenvalue",type="b");
    abline(v=nrFinal,col="blue");
    abline(h= 1,col="red")
    title(bquote(paste("Top ",.(nrFinal)," evs are ", .(round(100*frac[nrFinal],4)),
                       "% of total"))); #,cex=1.2);
    cat("Press Enter to continue...");
    readline();
    plot(ts, ur.h[,1], col=1, type="l", xlab="t", ylab="U", ylim=range(ur.h));
    for(i in 2:nrFinal)
      lines(ts, ur.h[,i], col=i);
    title(main = "showing the first few principle cpts")
  }

  zmat  <- lam.inv.sqrt  %*% t(ur.h) %*% (ys - mu.h);
  rv    <- list(des=des, t=ts, y=ys, nreps=nreps, ntps=ntps, mu=mu.h,
                sig=sig.h, esys=esys, ur.h=ur.h, zmat=zmat, nr=nrFinal
             );
  invisible(rv);

  
}  




########################################################
## function to emulate a list of points with projection
## ccs
########################################################
##
## emulate a list of points in the parameter space and
## return the back-projected functional results
##
## pointSet: matrix of points to be emulated
## (nrow -> nparams)
## (ncol -> number of points)
## thetas: vector of estimated hyperparams 
## (n -> nthetas -> (maybe) nparams + 2 )
## pca.decomp: pca decomposition generated by emu.pca or
## emu.pca.adapt
## 
## this is *substantially* faster than reconCurveAtPoint
## the setup and teardown for running the emulator through
## libEmu becomes a drag for Npoints > 10
##
## the results back from the C code are stored in the list
## emuResults, each column is the result from one of the z-emulators at each of
## the nEmupts
##
reconCurveAtList <- function(pointSet, thetas, pca.decomp, cov.fn.in=1, reg.order.in=1,
                             computeTCov=TRUE)
{
  nr <- dim(pca.decomp$zmat)[1]
  nparams <- dim(pca.decomp$des)[1]
  nmodelpts <- dim(pca.decomp$des)[2]
  # points = nemupts x nparams points
  nemupts <- dim(pointSet)[1]
  emuResult <- list(xpos = pointSet, mean = matrix(0, ncol=nr, nrow=nemupts), var=matrix(0, ncol=nr, nrow=nemupts))
  ntps <- length(pca.decomp$t)
  nthetas = length(thetas[1,])
  
  #print(emuResult$mean)
  
  for(i in 1:nr){
    ## oh my lord, the zmat index needs to change for the different models, or else, what happens?
    model <- list(xmodel=as.matrix(t(pca.decomp$des)), training=t(pca.decomp$zmat)[,i])

    ## 
    temp <- callEmulateAtList(model, as.numeric(thetas[i,]), as.numeric(pointSet),
                              nemupts, nmodelpts, nparams=nparams, nthetas=nthetas,
                              cov.fn=cov.fn.in, reg.order=reg.order.in)

    emuResult$mean[,i]<- temp$mean
    emuResult$var[,i] <- temp$var
  }


  #print(emuResult$mean)
  # now we need to recon the curve
  #
  # we read each result as a row of the emuResult$Mean matrix and
  # multiply it into the yMeanEmuRecon matrix, which is stored with each curve
  # in a single column
  yMeanEmuRecon <- matrix(0, ncol=nemupts, nrow=ntps)
  for(i in 1:nemupts){
    yMeanEmuRecon[,i] <- pca.decomp$mu + pca.decomp$ur %*% diag(sqrt(pca.decomp$esys$values[1:nr])) %*%
      emuResult$mean[i,]
  }


  ## this is a matrix of variance at a single location in the emulated
  ## set for each of the ntp cpts
  yVarEmuRecon <- matrix(0, ncol=nemupts, nrow=ntps)
  
  for(i in 1:ntps){
    for(j in 1:nemupts){
        # this is a bit confusing we're doing
      # yVar_i_j = Sum_k ( U_i_k * U_i_k * lambda_k * V_k_j)
      yVarEmuRecon[i,j] <- sum((pca.decomp$ur[i,]**2)*(pca.decomp$esys$values[1:nr])*emuResult$var[j,])
    }
  }
  

  finalResult <- list(xpos = pointSet, tvalues = pca.decomp$t, mean=yMeanEmuRecon, var = yVarEmuRecon)


  ## we may also want to compute the matrix of variances of each of the t cpts against
  ## each of the others, at each point in the emulated set
  ##
  ## we'll need to put a little bit of diagonal noise on this to stop the matrices from
  ## getting a small negative evalue, making them non pos def and blowing up the implaus
  if(computeTCov==TRUE){
    tCovList <- vector("list", nemupts)
    for(ptIndex in 1:nemupts){
      ## should be symmetric
      tVar <- matrix(0, nrow=ntps, ncol=ntps)
      for(i in 1:ntps){
        for(j in 1:ntps){
          tVar[i,j] <- sum(pca.decomp$ur[i,]*pca.decomp$ur[j,]*
            pca.decomp$esys$values[1:nr] * emuResult$var[ptIndex,])
        }
      }
      tCovList[[ptIndex]] <- tVar
    }
    finalResult$tCovList <- tCovList
  }
  invisible(finalResult)

}

