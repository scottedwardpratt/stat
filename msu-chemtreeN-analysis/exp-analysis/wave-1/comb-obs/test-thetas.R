## some functions to test what seems to be very variable output from the
## callEstimate routine
##
## we should compare the following for P.C estimated thetas and independently
## estimated data. Is the p.c somehow poorly conditioning things?
## 
## ccs, for a given observable index, we should:
## - histogram the distribution of thetas
## - do a par-coord plot of the thetas for each eval (against their index)
## - run the emulator at some fixed location and compare results (very important)

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

source("fnAnalysis.R")

## check the difference between the linear and constant regression results for a given obs index
plotLinearConst<- function (index=1){

nhtest <- 25
nobs <- 6

par(mfrow=c(2,1), mar=c(3,3,1,1))
load("thetas-comb-inde.dat")
plot(fn.list[[1]]$thetas[index,], col=index, xlab="index", ylab="theta-value", ylim=c(-2,5))
for(j in 1:nhtest){
  points(fn.list[[j]]$thetas[index,], col=index)
}
title("independent emulators")


load("fndata-thetas-test.dat")
plot(fn.list[[1]]$thetas[index,], col=index, xlab="index", ylab="theta-value", ylim=c(-2,5))
for(j in 1:nhtest){
  points(fn.list[[j]]$thetas[index,], col=index)
}
title("functional emulator")
}

## make a prediction at a fixed location for all the different sets of thetas
comparePredValues <- function(predLocn=c(10,50,0.05), thetafile="fndata-thetas-test.dat"){
  load(thetafile)
  ## now we have fn.list
  ntest <- length(fn.list)
  predlist <- vector("list", ntest)
  scale.vec <- attr(fn.list[[1]]$model.sample$des, "scaled:scale")
  cent.vec <- attr(fn.list[[1]]$model.sample$des, "scaled:center")
  predLocn.scaled <- t(as.matrix((predLocn - cent.vec) / scale.vec))

  meanMat <- matrix(0, nrow=ntest, ncol=6)
  varMat <- matrix(0, nrow=ntest, ncol=6)

  scale.vec <- attr(fn.list[[1]]$model.sample$y, "scaled:scale")
  cent.vec <- attr(fn.list[[1]]$model.sample$y, "scaled:center")

  
  for(i in 1:ntest){
    predlist[[i]] <- reconCurveAtList(predLocn.scaled, fn.list[[i]]$thetas.est, fn.list[[i]]$pca.decomp,
                                      cov.fn=fn.list[[i]]$cov.fn, reg.order=fn.list[[i]]$reg.order)
    meanMat[i,] <- predlist[[i]]$mean*scale.vec + cent.vec
    varMat[i,] <- predlist[[i]]$var*scale.vec + cent.vec
  }

  result <- list(mean=meanMat, var=varMat)
  result
}


## check the predictions
library("plotrix")
plotPredValues <- function(predvalues, index=1, exp.data=expData){
  plotCI(predvalues$mean[,index], uiw=1.96*sqrt(predvalues$var[,index]), ylab="y", xlab="run index")
  abline(h=exp.data$obsValue[index], col="red")
  abline(h=exp.data$obsValue[index] + exp.data$obsError[index], col="red", lty=2)
  abline(h=exp.data$obsValue[index] - exp.data$obsError[index], col="red", lty=2)
}

source("combAna.R")


# this is where the exp data was collected
locn.result <- c(10,50, 0.05)
# pick one of the training points and things work out fine
# this is just as well
locn.test <- designData[5,]
result.test <- c()
result.test$obsValue <- modelData[5,]
result.test$obsError <- sqrt(modelData[5,])
  
plist <- comparePredValues(predLocn=locn.result, thetafile="fndata-thetas-test.dat")
par(mfrow=c(2,3))
for(i in 1:6){
  plotPredValues(plist, index=i, exp.data=expData)
}
