## drop the grid to a csv file so we can
## look at it in paraview

dump.grid.csv <- function(data.grid, outname="grid-dump.csv"){
  npts <- dim(data.grid$grid)[1]

  xmin <- data.grid$xrange[1]
  ymin <- data.grid$yrange[1]
  zmin <- data.grid$zrange[1]

  xmax <- data.grid$xrange[2]
  ymax <- data.grid$yrange[2]
  zmax <- data.grid$zrange[2]

  # something was going wrong with using dx*npts to set the max
  # oh fuck because we need to add the offset of xmin first
  # ...
  if(is.null(xmax) || is.null(xmin) || is.null(ymax) || is.null(ymin) || is.null(zmax) || is.null(zmin)){
    stop("lattice file doesnt contain max min values, rebuild")
  }
  
  x.seq <- seq(xmin, xmax, length.out=npts)
  y.seq <- seq(ymin, ymax, length.out=npts)
  z.seq <- seq(zmin, zmax, length.out=npts)

  data.list.output <- expand.grid(x.seq, y.seq, z.seq)

  final.dat <- cbind(data.list.output, as.vector(data.grid$grid))

  colnames(final.dat) <- c("x", "y", "z", "implaus")

  write.csv(final.dat, file=outname, row.names=FALSE)
}
