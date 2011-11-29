##
## ccs, cec24@phy.duke,edu
## sep 2011
##i
## functions to setup and run the mc sampling process where we
## are trying to find probability distribution of parameters given
## a set of observations.
##
## to start we use the functions in fnAnalysis.R to generate
## a scaled sample, a pca decomposition of that sample, 
## and an estimated set of hyper-parameters for this model
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
##
## we then do metropolis sampling of density.observable(pt) * density.emulator(pt)
## and build up an empirical posterior probability distribution.
## we will need to define a function to compute density.observable(pt) and
## we sample the emulator with calls to fn.sampleEmu.density(fnData, pt)
##
## the emulator sampling requires a call to setup.MC.multi which allocates
## memory within libRbind so we don't have to invert too many covariance matrices
## when we're done we'll need to call destroy.MC.multi to free this memory

library(mvtnorm)
source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings

arch <- system("uname -s", intern=TRUE)
if(is.loaded("callEstimate") == FALSE){
  libNAME <- "~/local/lib/libRBIND"
  if(arch == "Linux"){
    libNAME <- paste(libNAME, ".so", sep="")
  } else if(arch =="Darwin"){
    libNAME <- paste(libNAME, ".dylib", sep="")
  } else {
    buffer <- paste("error: uname -s gives ", arch , " not supported", sep="")
    stop(buffer)
  }
  dyn.load(libNAME)
}

source("~/local/include/emu-analysis/fnAnalysis.R") # load the functional data helpers

## allocate the data we need to start sampling
setupSampling <- function(fn.data){
  nparams <- dim(fn.data$model.sample$des)[1]
  nydims <- fn.data$pca.decomp$nr

  if(fn.data$cov.fn==1){
    nthetas <- nparams + 2
  } else {
    nthetas <- 3
  }
   
  # setup the mc data structures
  cat("# setting up mc\n")
  setup.MC.multi(model=list(xmodel=t(fn.data$model.sample$des), training=t(fn.data$pca.decomp$zmat)),
                 thetas=fn.data$thetas.est,
                 nmodelpoints=fn.data$nruns,
                 nparams=nparams, 
                 nydims=nydims,
                 nthetas=nthetas,
                 cov.fn=fn.data$cov.fn,
                 reg.order=fn.data$reg.order)

}

## @param x0 initial vector location in the design space
## @param z0 initial vector location in the measurement space
## @param data.density.fn a function which returns the log(density) of the observational data at a location in z
## 
sample.integrand <- function(fn.data, exp.data, nsteps=10000, x0, z0, xranges, zranges, data.density.fn){
  setupSampling(fn.data)

  naccept <- 0
  ndims.x <- length(x0)
  ndims.z <- length(z0)
  ndims <- ndims.x + ndims.z

  x.prev <- x0
  z.prev <- z0
  delta.x <- xranges[,2] - xranges[,1]
  delta.z <- zranges[,2] - zranges[,1]

  samples <- matrix(0, nrow=nsteps, ncol=ndims)

  ## a guess for the stepping scale 
  alpha.x <- rep(0.01, ndims.x)
  alpha.z <- rep(0.01, ndims.z)



  ## main mc sampling loop
  for(i in 1:nsteps){

    if(i %% (nsteps/10) == 0){
      pacc <- naccept / i
      cat("# sampled: ", i, " paccept: ", pacc, "\n")
    }

    ## generate a proposal step in the parameter space
    x.star <- x.prev + alpha.x*delta.x*runif(ndims.x, min=-1, max=1)
    ## check bounds for x.star
    for(index in 1:ndims.x){
      if(x.star[index] > xranges[index, 2] || x.star[index] < xranges[index, 1]){
        x.star[index] <- xranges[index,1] + (x.star[index] - xranges[index,1])%%delta.x[index]
      }
    }

    ## generate a normally distributed proposal step in the observable space
    #z.star <- rnorm(ndims.z, z.prev, alpha.z)
    z.star <- z.prev + alpha.z*delta.z*runif(ndims.z, min=-1, max=1)
    ## check bounds for x.star
    for(indez in 1:ndims.z){
      if(z.star[indez] > zranges[indez, 2] || z.star[indez] < zranges[indez, 1]){
        z.star[indez] <- zranges[indez,1] + (z.star[indez] - zranges[indez,1])%%delta.z[indez]
      }
    }


    ## compute the acceptance ratio for this step
    ## sample the emulator and the data.density.fn
    ## for the numerator and denominators
    log.r.num <- fn.sampleEmu.density(fn.data, x.star, z.star)
    + data.density.fn(exp.data, z.star)

    log.r.denom <- fn.sampleEmu.density(fn.data, x.prev, z.prev)
    + data.density.fn(exp.data, z.prev)

    ## now include the density due to the prior over the design space
    for(index in 1:(ndims.x-1)){
      ## sample a uniform prior on the design space
      log.r.num <- log.r.num + dunif(x.star[index], min=xranges[index,1], max=xranges[index,2], log=TRUE)
      log.r.denom <- log.r.denom + dunif(x.prev[index], min=xranges[index,1], max=xranges[index,2], log=TRUE)
    }

    log.r <- log.r.num - log.r.denom

    
    if(log(runif(1)) < log.r){
      x.new <- x.star
      z.new <- z.star
      naccept <- naccept + 1
    } else {
      x.new <- x.prev
      z.new <- z.prev
    }

    samples[i,] <- c(x.new, z.new)
    x.prev <- x.new
    z.prev <- z.new
    
  }
  ## clear up
  nydims.emu <- fn.data$pca.decomp$nr
  destroy.MC.multi(nydims.emu)

  cat("# p(accept): ", naccept/nsteps, "\n")
  samples
}



rescale.samps <- function(samples, fn.data){
  nparams <- dim(fn.data$model.sample$des)[1]
  nbins <- fn.data$nbins
  nsamps <- dim(samples)[1]
  scale.vec <- c(attr(fn.data$model.sample$des, "scaled:scale"), attr(fn.data$model.sample$y, "scaled:scale"))
  center.vec <- c(attr(fn.data$model.sample$des, "scaled:center"), attr(fn.data$model.sample$y, "scaled:center"))

  samples.rescaled <- matrix(0, nrow=nsamps, ncol=(nparams+nbins))

  for(i in 1:nsamps)
    samples.rescaled[i,] <- samples[i,] * scale.vec + center.vec

  samples.rescaled
}



#sample.integrand <- function(fn.data, exp.data, nsteps=10000, x0, z0, xranges, zranges, data.density.fn){

## assume that we're on the unit interval here
generate.mc.samples <- function(fn.data, exp.data, data.density.fn, nprocs=4, nstep.mc=500000){
  library("multicore")
  samp.par <- vector("list", nprocs)

  ic.x <- matrix(runif(nparams*nprocs, 0, 1), nrow=nprocs, ncol=nparams)
  ic.z <- matrix(runif(nbins*nprocs, -pi/2, pi/2), nrow=nprocs, ncol=nbins)

  xr.in <- matrix(c(0,1), nrow=nparams, ncol=2, byrow=TRUE)
  ## we'll assume that all the exp distributions are normal and so on +- infinity
  ## we can sample this by taking our uniform range as +- pi/2 and then doing an atan when
  ## we actually compute the density
  zr.in <- matrix(c(-pi/2,pi/2), nrow=nbins, ncol=2, byrow=TRUE)

  ## run the sampling nprocs at a time
  for(i in 1:nprocs){
    samp.par[[i]] <- parallel(sample.integrand(fnData, exp.data, nsteps=nstep.mc,
                                                    x0=ic.x[i,], z0=ic.z[i,],
                                                    xranges=xr.in,
                                                    zranges=zr.in,
                                                    data.density.fn),
                              mc.set.seed=TRUE)
    
  }
  ## collate the parallel result
  res <- collect(samp.par)
  res.burnt <- vector("list", nprocs)
  
  ## currently we'll throw away the first fifth as  burn in
  nburn <- nstep.mc / 5
  for(i in 1:nprocs){
    res.burnt[[i]] <- res[[i]][nburn:nstep.mc,]
  }
  stats <- mc.stats(res.burnt)

  samples.rescaled <- vector("list", nprocs)
  for(i in 1:nprocs){
    samples.rescaled[[i]] <- rescale.samps(res.burnt[[i]], fnData)
  }

  samples.final <- NULL
  for(i in 1:nprocs){
    samples.final <- rbind(samples.final, samples.rescaled[[i]])
  }

  result <- list(samples=samples.final, stats=stats)
  invisible(result)
}



## compute the within and between run stats for
## some list of j sets of n observations samples.par
## where psi_ij is the ith oversation of the j'th process
##
## B = n / (j-1) sum_j (psi_j - psi_tot)^2
## the between process variance
##
## where psi_j = 1/n sum_i psi_ij
## psi_tot = 1/J sum_j psi_j
##
## W = 1/J sum_j S_j^2
## the within process variance
## 
## S_j^2 = 1/(n-1) Sum_i (psi_ij - psi_j) **2
mc.stats <- function(samples.par, obs.index=1){
  ## how many parallel sets of runs we have
  nsets <- length(samples.par)
  ## how many instances we have for each run
  npts <- dim(samples.par[[1]])[1]

  mean.vec <- rep(NA, nsets)
  for(i in 1:nsets){
    # we'll just average EVERYTHING in the sample
    mean.vec[i] <- mean(samples.par[[i]][,obs.index])
  }
  mean.tot <- mean(mean.vec)
  
  ## this is the variance of the means of each process
  B.var <- 0
  for(i in 1:nsets){
    B.var <- B.var + (mean.vec[i] - mean.tot)**2
  }
  B.var <- B.var *(npts / (nsets-1))

  ## compute the within process variance for each proc
  ## this is just the mean of the variance of each set of samples
  S.vec <- rep(NA, nsets)
  for(i in 1:nsets){
    S.vec[i] <- var(samples.par[[i]][,obs.index])
  }
  ## 
  W.var <- mean(S.vec)

  ## the estimated posterior variance (over est)
  ## var( psi | y )
  Var.post <- (npts - 1)/(npts) * W.var + 1/npts * B.var

  R.scale <-  sqrt(Var.post / W.var)

  result <- list(B.var=B.var, W.var=W.var, Var.post.est=Var.post, R.scale = R.scale)
  invisible(result)
}


