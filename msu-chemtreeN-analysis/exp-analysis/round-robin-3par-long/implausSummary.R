##
## ccs, cec24@phy.duke.edu
## 25.01.2012
##
## compute some simple integrals of the grids and see if we can use that to compare

l2.dist.grid <- function(grid1, grid2){
  if(grid1$dx != grid2$dx || grid1$dy != grid2$dy || grid1$dz != grid2$dz){
    dvol <- (mean(grid1$dx,grid2$dx)*mean(grid1$dy,grid2$dy)*mean(grid1$dz*grid2$dz))
  } else {
    dvol <- grid1$dx*grid1$dy*grid1$dz ## add up the dx bits
  }
  #cat("# dvol: ", dvol, "\n")
  
  npts <- dim(grid1$grid)
  product <- grid1$grid * grid2$grid

  vol <- (grid1$xrange[2]-grid1$xrange[1])*(grid1$yrange[2]-grid1$yrange[1])*(grid1$zrange[2]-grid1$zrange[1])

  norm <- sum(product*dvol) * vol

  norm
}

l2.dist.grid.inv <- function(grid1, grid2){
  if(grid1$dx != grid2$dx || grid1$dy != grid2$dy || grid1$dz != grid2$dz){
    dvol <- (mean(grid1$dx,grid2$dx)*mean(grid1$dy,grid2$dy)*mean(grid1$dz*grid2$dz))
  } else {
    dvol <- grid1$dx*grid1$dy*grid1$dz ## add up the dx bits
  }
  #cat("# dvol: ", dvol, "\n")
  
  npts <- dim(grid1$grid)
  product <- (1/grid1$grid) * (1/grid2$grid)

  norm <- sum(product*dvol)

  norm
}


## symmetric normalized metric? basically computes the projection of grid1 onto grid 2
l2.normed.grid <- function(grid1, grid2){
  g1.l2 <- l2.dist.grid(grid1, grid1)
  g2.l2 <- l2.dist.grid(grid2, grid2)

  g12 <- l2.dist.grid(grid1, grid2)

  g12 / (sqrt(g1.l2) * sqrt(g2.l2))
}

l2.normed.grid.inv <- function(grid1, grid2){
  g1.l2 <- l2.dist.grid.log.inv(grid1, grid1)
  g2.l2 <- l2.dist.grid.log.inv(grid2, grid2)

  g12 <- l2.dist.grid.log.inv(grid1, grid2)

  g12 / (sqrt(g1.l2) * sqrt(g2.l2))
}



## assymm metric, facundo's suggestion, use the trained on matrix as the normalizing case
l2.normed.grid.train <- function(grid.comp, grid.train){
  #gcomp.l2 <- l2.dist.grid(grid.comp, grid.comp)
  gtrain.l2 <- l2.dist.grid(grid.train, grid.train)

  g12 <- l2.dist.grid(grid.comp, grid.train)

  ## return the ij distance normalized by grid 1
  g12 / gtrain.l2
}


## subtract the two distributions
int.sub.grids <- function(grid1, grid2){
  if(grid1$dx != grid2$dx || grid1$dy != grid2$dy || grid1$dz != grid2$dz){
    dvol <- (mean(grid1$dx,grid2$dx)*mean(grid1$dy,grid2$dy)*mean(grid1$dz*grid2$dz))
  } else {
    dvol <- grid1$dx*grid1$dy*grid1$dz ## add up the dx bits
  }
  #cat("# dvol: ", dvol, "\n")
  
  npts <- dim(grid1$grid)
  diff <- (grid1$grid - grid2$grid)^2

  norm <- sum(diff*dvol)
}

int.sum.grids <- function(grid1, grid2){
  if(grid1$dx != grid2$dx || grid1$dy != grid2$dy || grid1$dz != grid2$dz){
    dvol <- (mean(grid1$dx,grid2$dx)*mean(grid1$dy,grid2$dy)*mean(grid1$dz*grid2$dz))
  } else {
    dvol <- grid1$dx*grid1$dy*grid1$dz ## add up the dx bits
  }
  #cat("# dvol: ", dvol, "\n")
  
  npts <- dim(grid1$grid)
  sum.grid <- (grid1$grid + grid2$grid)^2

  norm <- sum(sum.grid*dvol)
}

## compute the fraction of the grid with values less than thresh
int.grid.thresh <- function(grid, thresh){
    dvol <- grid$dx*grid$dy*grid$dz ## add up the dx bits

    tvals <-  sapply(grid$grid, function(x)(if(x<thresh){x}else{0.0}))
    integral.part <- sum(tvals*dvol)

    int.full <- sum(dvol*grid$grid)

    integral.part / int.full
}

norm.sub.grids <- function(grid.comp, grid.train){
  diff <- int.sub.grids(grid.comp, grid.train) ## do we want to normalize this in any way?
  norm.train <- l2.dist.grid(grid.train, grid.train)
  norm.comp <- l2.dist.grid(grid.comp, grid.comp)

  diff / norm.train
}

norm.sum.diff.grids <- function(grid1, grid2){
  diff <- int.sub.grids(grid1, grid2)
  sum <- int.sum.grids(grid1, grid2)

  sqrt(diff / sum)
}


compare.round.robin <- function(outfile="./l2-comparison-matrix.dat"){
  subfolders <- c("mw1", "mw2", "mw3", "mw6")
  topfolders <- c("mw2", "mw3", "mw6")

  ## mw2Grids <- vector("list", 4)
  ## mw3Grids <- vector("list", 4)
  ## mw6Grids <- vector("list", 4)

  ngrids <- length(subfolders)*length(topfolders)
  
  gridList <- vector("list", ngrids)
  trainingGrids <- vector("list", length(topfolders))

  normMat <- matrix(0, nrow=ngrids, ncol=ngrids)
  normMat.inv <- matrix(0, nrow=ngrids, ncol=ngrids)
  sum.diff.Mat <- matrix(0, nrow=ngrids, ncol=ngrids)
  ## compute the l2 norm against the trained emulator only
  ## 
  normMat.train <- matrix(0, nrow=length(topfolders), ncol=ngrids)
  ## compute integral of (f - g)^2 over the whole grid where we always compare against the training data
  ## 
  diffMat.train <- matrix(0, nrow=length(topfolders), ncol=ngrids)

  count <- 1
  count.train <-  1
  thresh.seq <- seq(0,30,length.out=60)
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
      #gridList[[count]]$thresh <- sapply(thresh.seq, function(x)(int.grid.thresh(implausGrid, x)))
      if(mwString == mwStringCompare){
        trainingGrids[[count.train]] <- implausGrid
        trainingGrids[[count.train]]$name <- grid.Name
        count.train <- count.train + 1
      }

      count <- count + 1
    }
  }

  
  names <- sapply(gridList, function(x)(x$name))

  colnames(normMat) <- names
  rownames(normMat) <- names

  colnames(normMat.inv) <- names
  rownames(normMat.inv) <- names

  ## thresh.vals <- sapply(gridList, function(x)(x$thresh))
  ## thresh.final <- cbind(thresh.seq, thresh.vals)
  ## colnames(thresh.final) <- c("t", names)

  ## #browser()
  ## write.table(thresh.final, file="thresh-matrix.dat", row.names=FALSE)

  
  
  #browser()
  
  for(i in 1:ngrids){
    for(j in 1:ngrids){
      cat("# ", i, j, "\n")
      normMat[i,j] <- l2.normed.grid(gridList[[i]], gridList[[j]])
      #normMat.inv[i,j] <- l2.normed.grid.inv(gridList[[i]], gridList[[j]])
      sum.diff.Mat[i,j] <- norm.sum.diff.grids(gridList[[i]], gridList[[j]])
    }
  }

  rownames(normMat.train) <- sapply(trainingGrids, function(x)(x$name))
  colnames(normMat.train) <- names

  rownames(diffMat.train) <- sapply(trainingGrids, function(x)(x$name))
  colnames(diffMat.train) <- names
  
  for(i in 1:length(topfolders)){
    for(j in 1:ngrids){
      normMat.train[i,j] <- l2.normed.grid.train(gridList[[j]], trainingGrids[[i]])
      diffMat.train[i,j] <- norm.sub.grids(gridList[[j]], trainingGrids[[i]])
    }
  }



  write.table(normMat, file=outfile)
  write.table(normMat.inv, file="l2-norm-inv.dat")
  write.table(normMat.train, file="l2-modified-training-compare.dat")
  write.table(diffMat.train, file="sub-training-compare.dat")
  write.table(sum.diff.Mat, file="diff-sum-grid-compare.dat")
  
  invisible(normMat)
}


plot.thresh <- function(infile="thresh-matrix.dat"){
  library(sp)
  library(lattice)
  tm <- read.table(infile)
  xyplot( obs.mw3.train.mw3 + obs.mw2.train.mw2 + obs.mw6.train.mw6 ~ t, tm, auto.key=list(title="obs", space="right", cex=1.0))
}
 
