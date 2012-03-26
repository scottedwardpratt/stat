## drop the grid to a csv file so we can
## look at it in paraview

dump.grid.csv <- function(data.grid, outname="grid-dump.csv", xmin=0, ymin=0, zmin=0){
  npts <- dim(data.grid$grid)[1]

  x.seq <- seq(xmin, data.grid$dx*npts, length.out=npts)
  y.seq <- seq(ymin, data.grid$dy*npts, length.out=npts)
  z.seq <- seq(zmin, data.grid$dz*npts, length.out=npts)

  data.list.output <- expand.grid(x.seq, y.seq, z.seq)

  browser()
  
  final.dat <- cbind(data.list.output, as.vector(data.grid$grid))

  colnames(final.dat) <- c("x", "y", "z", "implaus")

  write.csv(final.dat, file=outname)
}
