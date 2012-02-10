## redump grid files to csv from .dat files
## ccs cec24@phy.duke.edu
## 10.02.2012
##
## the .dat files are correct for all the grids (yay) but the
## grid.dump method was using the wrong ranges because i failed to add the minimum values of the
## variables to offset the range to compute the max.
## it was:  xmax = npts*dx
## it is now: xmax = npts*dx + xmin
## 
## this files zoops through the folders and redumps the grids
source("dump-grid-csv.R")

roundFolders <- c("./mw2-round", "./mw3-round", "./mw6-round")
roundStubs <- c("mw2", "mw3", "mw6")
subFolders <- c("mw1", "mw2", "mw3", "mw6")

for( i in 1:3){
  round <- roundFolders[i]
  roundStub <- roundStubs[i]
  for(j in 1:4){
    sub <- subFolders[j]
    buffer.dat <- paste(round,"/", sub, "/implaus-grid-observations-", sub, "-training-", roundStub, ".dat", sep="")
    buffer.out <- paste(round, "/", sub, "/implaus-grid-observations-", sub, "-training-", roundStub, "-dump.csv", sep="")    
    cat("# processing: ", buffer.dat, "\n")
    cat("# outputting: ", buffer.out, "\n")
    load(buffer.dat)
    dump.grid.csv(implausGrid, outname=buffer.out)
  }
}
