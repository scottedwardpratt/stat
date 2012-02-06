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


  norm <- sum(product*dvol)

  norm
}

l2.dist.grid.log.inv <- function(grid1, grid2){
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

l2.normed.grid.log.inv <- function(grid1, grid2){
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

norm.sub.grids <- function(grid.comp, grid.train){
  diff <- int.sub.grids(grid.comp, grid.train) ## do we want to normalize this in any way?
  norm.train <- l2.dist.grid(grid.train, grid.train)
  norm.comp <- l2.dist.grid(grid.comp, grid.comp)

  diff / norm.train
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
  normMat.log.inv <- matrix(0, nrow=ngrids, ncol=ngrids)
  ## compute the l2 norm against the trained emulator only
  ## 
  normMat.train <- matrix(0, nrow=length(topfolders), ncol=ngrids)
  ## compute integral of (f - g)^2 over the whole grid where we always compare against the training data
  ## 
  diffMat.train <- matrix(0, nrow=length(topfolders), ncol=ngrids)

  count <- 1
  count.train <-  1
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

  colnames(normMat.log.inv) <- names
  rownames(normMat.log.inv) <- names

  
  #browser()
  
  for(i in 1:ngrids){
    for(j in 1:ngrids){
      cat("# ", i, j, "\n")
      normMat[i,j] <- l2.normed.grid(gridList[[i]], gridList[[j]])
      normMat.log.inv[i,j] <- l2.normed.grid.log.inv(gridList[[i]], gridList[[j]])
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
  write.table(normMat.log.inv, file="l2-norm-log-inv.dat")
  write.table(normMat.train, file="l2-modified-training-compare.dat")
  write.table(diffMat.train, file="sub-training-compare.dat")
  
  invisible(normMat)
}

