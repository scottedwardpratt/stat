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

source("fnMcSampling.R")

nsamples <- 200

## load the model data
##
## first the luminosity data
lumOutputFile  <- "./wave-1-long/lum_fun_outputs_2.dat"
modelDataLum <- as.matrix(read.table(lumOutputFile))
nbinsLum <- dim(modelDataLum)[2]
nruns <- dim(modelDataLum)[1]
##
## now the metallicity data
metOutputFile <- "./wave-1-long/metallicity_MV_outputs_2.dat"
modelDataMet <- as.matrix(read.table(metOutputFile))
nbinsMet <- dim(modelDataMet)[2]

if(nruns != dim(modelDataMet)[1]){
  stop("nruns modelDataMet doesn't match modelDataLum")
}

## redefine nruns
nruns <- nsamples


nbins <-  nbinsMet + nbinsLum 
#nbins <- nbinsLum
modelData.big <- cbind(modelDataLum, modelDataMet)

modelData <- modelData.big

## load the design
designFile <- "./design/design_2_sorted.dat"
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designData.big <- as.matrix(read.table(designFile, col.names=desNames))

designData <- designData.big


## load the experimental data
##
## lum, need to use the _2 observations with the long data set
expDataFileLum <- "./lum_fun_observed_2.dat"
expDataLum <- as.matrix(read.table(expDataFileLum))

expDataFileMet <- "./metallicity_MV_observed_2.dat"
expDataMet <- as.matrix(read.table(expDataFileMet))

## the lum data errors are "poisson" ~ sqrt(\lambda), if N_obs was very big we could approximate it as a normal
## distribution with mean \lambda and sd \sqrt{\lambda}
## 
## the met data are 95% confidence bounds, 95% of the data fall within mu +- 2 \sigma
## so sigma = upper-lower / 2 
## where upper = mean + conf, lower = mean - conf
## so sigma = ((mean + conf ) - (mean - conf))  / 2 = conf
expData <- list(obsValue=c(expDataLum[2,], expDataMet[1,]),
                  obsError=c(expDataLum[3,], expDataMet[2,]))

rebuild <- 0
buffer <- paste("functional-data-lum-", nsamples, ".dat", sep="")
if(rebuild == 1 || file.exists(buffer) == FALSE){
  ##
  ## generate a functional sample from the vars in global sope
  fnData <- fn.sample.gen(cov.fn=1, reg.order=1)
  ## now do the pca decomp
  fnData <- fn.pca.gen(fnData, cutOff=0.99)
  ## estimate the thetas
  fnData <- fn.estimate(fnData)

  save(fnData, file=buffer)
} else {
  load(buffer)
}

## we need to define a function that will give a log-density for our experimental data (exp.data)
## given a particular proposed value (z.prop)
## 
## for this test lets pretend that all our experimental data has a
## normal distribution
##
## BUT!
## the emulator is still scaled to the unit cube during the mc sampling, so we need
## to scale the exp.data correctly
##
## AND!
## the support of dnorm is +- infinity, we need to draw z.prop over this range, we
## draw z.prop' on +- pi/2 and then atan it to get it onto +- infty
exp.density.fn <- function(exp.data.scaled, z.prop){
  nzdim <- length(z.prop)
  density <- rep(0, nzdim)
  ## we sample z.prop on +- pi/2 and then unwrap this to -+ infty using tan(z.prop) 
  z.prop <- tan(z.prop)
  for(i in 1:nzdim){
    density[i] <- dnorm(z.prop[i], mean=exp.data.scaled$obsValue[i], sd=exp.data.scaled$obsError[i], log=TRUE)
  }
  sum(density)
}

## scale the observables so that we can use them in the mc sampling
scale.obs <- function(fn.data, exp.data){
  scale.vec.obs <- attr(fn.data$model.sample$y, "scaled:scale")
  center.vec.obs <- attr(fn.data$model.sample$y, "scaled:center")
  result <- exp.data

  result$obsValue <- (result$obsValue - center.vec.obs)/(scale.vec.obs)
  ## not sure is the error a sd or a var? it's a var...
  result$obsError <- (result$obsError)/(scale.vec.obs)

  result
}

expData.scaled <- scale.obs(fnData, expData)

samples.final <- generate.mc.samples(fnData, expData.scaled, exp.density.fn, nstep.mc=5000000)

## dump the mc.samples for later
save(samples.final, file="functional-samples-200.dat")
#load("functional-samples.dat")

plotSamples <- function(samples){
  library(MASS)
  ndims.x <- 3
  ## one way to look at the data, not very exciting
  #pairs(samples[,1:ndims.x], cex=0.1, labels=desNames)

  smoothFn <- function(x,y){
    density <- kde2d(x, y)
    probs <-c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995)
    qs <- c(0, quantile(density$z, probs))
    nlevels <- length(qs)
    breaks <- c(qs, max(density$z))
    image(density, col=cm.colors(nlevels), breaks=breaks,
          add=TRUE)
    contour(density, add=TRUE, levels=qs, labels=c(0, probs))
  }
  
  ## make 2d kernel smoother estimates for each pair
  pairs(samples[,1:ndims.x], panel=smoothFn, labels=desNames)

  
}

plotSamples(samples.final$samples)
