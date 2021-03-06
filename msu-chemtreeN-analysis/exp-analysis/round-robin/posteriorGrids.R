##
## ccs, cec24@phy.duke.edu
## use implausibility grids to compute bayesian evidence for the combination of
## model (training data) and observable
##
## evidence(M) = Int P(d | theta, M) * P(theta, M) dtheta
##
## where theta -> parameters of model
## d -> model observations
## P(theta, M) -> prior (we'll take these to be uniform)
## P(d |theta, M) -> posterior, we can use Exp(-I/2) since we're actually reporting I^2 in the code and then
## we have somethign like a chi^2 dist
## 
##
compute.evidence <- function(grid){
  dvol <- grid$dx * grid$dy * grid$dz

  ## take the prior as uniform
  priorValueX <- 1.0/(grid$xrange[2]-grid$xrange[1])
  priorValueY <- 1.0/(grid$yrange[2]-grid$yrange[1])
  priorValueZ <- 1.0/(grid$zrange[2]-grid$zrange[1])
  
  
  npts <- dim(grid$grid)

  posterior <- exp(-grid$grid/2.0)

  ## note that the true volume of the space and the prior exactly cancel out
  evidence <- sum(posterior * dvol)
}

compute.evidence.pt <- function(grid, pt){
  xbin <- round((pt[1] - grid$xrange[1]) / grid$dx)
  ybin <- round((pt[2] - grid$yrange[1]) / grid$dy)
  zbin <- round((pt[3] - grid$zrange[1]) / grid$dz)

  gridValue <- grid$grid[xbin, ybin, zbin]

  posterior <- exp(-gridValue / 2.0) 


  posterior
}
  


library(lattice)
library(sp)
plot.posterior <- function(grid, nslices=9){
  trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()

  posterior <- exp(-grid$grid / 2.0)

  npts <- dim(grid$grid)[1]

  x <- seq(from=grid$xrange[1], to=grid$xrange[2],length.out=npts)
  y <- seq(from=grid$yrange[1], to=grid$yrange[2],length.out=npts)
  
  des.grid <- expand.grid(Zr=x, Fe=y)

  dzSlice <- (grid$zrange[2] - grid$zrange[1]) / nslices
  dzIndex <- floor(npts / nslices)
  zind <- 1

  sliceNames <- c()
  
  for(index in 1:nslices){

    zval <- index*dzSlice # cts z val
    zind <- index*dzIndex
    buffer <- paste("Fb", round(zval,2), sep="")
    sliceNames <- c(sliceNames, buffer)
    cat("buffer: ", buffer, "\n")
    des.grid[[buffer]] <- as.vector(posterior[,,zind])
  }

  namesPlus <- paste(sliceNames, collapse=" + ", sep="")
  cat("namesPlus: ", namesPlus, "\n")
  eqn <- as.formula(paste(namesPlus, "~ Zr * Fe", sep=""))

  p1 <- levelplot(eqn, des.grid, main="posterior density", scales=list(cex=0.8), layout=c(3,3))

  p1
}

plot.all.posts <- function(){
  subfolders <- c("mw1", "mw2", "mw3", "mw6")
  topfolders <- c("mw2", "mw3", "mw6")
  for(i in 1:length(topfolders)){
    for(j in 1:length(subfolders)){
      mwString <- topfolders[i]
      mwStringCompare <- subfolders[j]
      buffer <- paste("./", mwString, "-round/", mwStringCompare , "/implaus-grid-observations-", mwStringCompare, "-training-", mwString, ".dat", sep="")
      cat("# ", buffer, "\n")
      load(buffer)
      p <- plot.posterior(implausGrid)
      buffer <- paste("./", mwString, "-round/", mwStringCompare , "/images/posterior-", mwStringCompare, "-training-", mwString, ".pdf", sep="")
      pdf(buffer)
      trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()
      print(p)
      dev.off()
    }
  }
}


compare.evidence.at.point <- function(point){
  subfolders <- c("mw1", "mw2", "mw3", "mw6")
  topfolders <- c("mw2", "mw3", "mw6")

  ngrids <- length(subfolders)*length(topfolders)
                         
  gridList <- vector("list", ngrids)

  evidMat <- matrix(0, nrow=ngrids, ncol=ngrids)

  count <- 1
  for(i in 1:length(topfolders)){
    for(j in 1:length(subfolders)){
      mwString <- topfolders[i]
      mwStringCompare <- subfolders[j]
      buffer <- paste("./", mwString, "-round/", mwStringCompare , "/implaus-grid-observations-", mwStringCompare, "-training-", mwString, ".dat", sep="")
      cat("# ", buffer, "\n")
      load(buffer)

      grid.Name <- paste("obs-",mwStringCompare, "-train-", mwString, sep="")
      
      gridList[[count]] <- implausGrid
      gridList[[count]]$name <- grid.Name
      count <- count + 1
    }
  }

  names <- sapply(gridList, function(x)(x$name))
  colnames(evidMat) <- names
  rownames(evidMat) <- names

  
  for(i in 1:ngrids){
    for(j in 1:ngrids){
      e1 <- compute.evidence.pt(gridList[[i]], point)
      e2 <- compute.evidence.pt(gridList[[j]], point)

      cat("# ", i, j, e1, e2, "\n")
      if(i != j){
        evidMat[i,j] <- log10(e1/e2)
      } else {
        evidMat[i,j] <- log10(e1/e1)
      }
    }
  }

  evidMat
}
  
  


compare.evidence <- function(){
  subfolders <- c("mw1", "mw2", "mw3", "mw6")
  topfolders <- c("mw2", "mw3", "mw6")

  ngrids <- length(subfolders)*length(topfolders)
                         
  gridList <- vector("list", ngrids)

  evidMat <- matrix(0, nrow=ngrids, ncol=ngrids)

  count <- 1
  for(i in 1:length(topfolders)){
    for(j in 1:length(subfolders)){
      mwString <- topfolders[i]
      mwStringCompare <- subfolders[j]
      buffer <- paste("./", mwString, "-round/", mwStringCompare , "/implaus-grid-observations-", mwStringCompare, "-training-", mwString, ".dat", sep="")
      cat("# ", buffer, "\n")
      load(buffer)

      grid.Name <- paste("obs-",mwStringCompare, "-train-", mwString, sep="")
      
      gridList[[count]] <- implausGrid
      gridList[[count]]$name <- grid.Name
      count <- count + 1
    }
  }

  names <- sapply(gridList, function(x)(x$name))
  colnames(evidMat) <- names
  rownames(evidMat) <- names

  
  for(i in 1:ngrids){
    for(j in 1:ngrids){
      e1 <- compute.evidence(gridList[[i]])
      e2 <- compute.evidence(gridList[[j]])

      cat("# ", i, j, e1, e2, "\n")
      if(i != j){
        evidMat[i,j] <- log10(e1/e2)
      } else {
        evidMat[i,j] <- log10(e1/e1)
      }
    }
  }

  evidMat
}

## on source this will make all the posterior plots
##
plot.all.posts()
## and save the evidence matrix
emat <- compare.evidence()
## save as a text table, with row and column labels hopefully
write.table(emat, file="evidence-matrix.dat")

## the true location
pt <- c(10, 50, 0.05)

emat.pt <- compare.evidence.at.point(pt)
write.table(emat.pt, file="evidence-matrix-at-pt.dat")
