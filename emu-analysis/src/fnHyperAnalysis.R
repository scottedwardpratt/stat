library(lhs)

## each location in the plane is the result of averaging
## nsamples locations drawn by LHS samples in the orthogonal directions
## to the cube we're constructing
## 
fn.emulate.hyper.slice <- function(fn.data, dimA, dimB, stepDim, stepValue, nemupts=32, nsamples=64){
  des.scaled <- t(fn.data$model.sample$des)
  ## now create the correct ranges for sampling the slice
  rangeA <- c(min(des.scaled[,dimA]), max(des.scaled[,dimA]))
  rangeB <- c(min(des.scaled[,dimB]), max(des.scaled[,dimB]))
  stepSizeA <- (rangeA[2] - rangeA[1]) / (nemupts)
  stepSizeB <- (rangeA[2] - rangeA[1]) / (nemupts)

  noSample.vec <- c(dimA, dimB, stepDim) ## these are the dirns not to sample in
  sample.vec <- seq(1,nparams)[-noSample.vec]
  
  min.des <- apply(fn.data$model.sample$des, 1, min)[-noSample.vec]
  max.des <- apply(fn.data$model.sample$des, 1, max)[-noSample.vec]
  
  ## the points to sample at
  nsample.dirs <- nparams - 3
  pointList <- matrix(0, nrow=nsamples, ncol=nparams)

  meanList <- vector("list", nemupts**2)
  varList <- vector("list", nemupts**2)
  tCovList <- vector("list", nemupts**2)
  
  for(i in 1:nemupts){
    for(j in 1:nemupts){
      pointVec.A <- rep(rangeA[1]+i*stepSizeA, nsamples)
      pointVec.B <- rep(rangeB[1]+j*stepSizeB, nsamples)
      stepVec <- rep(stepValue, nsamples)

      ## now we generate points in the remaining dimensions
      samples <- maximinLHS(nsamples, nsample.dirs)*(max.des-min.des) + min.des

      pointList[,dimA] <- pointVec.A
      pointList[,dimB] <- pointVec.B
      pointList[,stepDim] <- stepVec
      pointList[,sample.vec] <- samples

      ## clist = list(xpos, tvalues, mean, var)
      ## where xpos == pointList, tvalues == pca.decomp$t
      clist <- reconCurveAtList(pointList, fn.data$thetas.est, fn.data$pca.decomp,
                                cov.fn=fn.data$cov.fn, reg.order=fn.data$reg.order)

      meanList[[j+nemupts*(i-1)]] <- clist$mean
      varList[[j+nemupts*(i-1)]] <- clist$var
      tCovList[[j+nemupts*(i-1)]] <- clist$tCovList
    }
  }

  result <- list(
                 mean=meanList, 
                 var=varList,
                 tCovList=tCovList,
                 reg.order=fn.data$reg.order, cov.fn=fn.data$cov.fn,
                 dimA=dimA, dimB=dimB, stepDim=stepDim, stepValue=stepValue,
                 npts=nemupts)
  
}

fn.reduce.hyper.slice <-  function(slice){
  npts <- slice$npts
  mean <- matrix(0, nrow=npts**2, ncol=nbins)
  var <- matrix(0, nrow=npts**2, ncol=nbins)
  
  for(i in 1:npts){
    for(j in 1:npts){
      mean <- apply(slice$mean[[j+npts*(i-1)]]$mean, 1, mean)
      var <- apply(slice$mean[[j+npts*(i-1)]]$var, 1, mean)
    }
  }
  result <- list(mean=mean, var=var, tCovList=slice$tCovList, reg.order=slice$reg.order,
                 cov.fn=slice$cov.fn, npts=npts)
}

