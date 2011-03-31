source("gen-samp.R") # load the sample set
source("emu-pca.R") # load the pca fns

# note, need to change this to .so or whatever for deployment
dyn.load("~/local/lib/libRBIND.dylib") #
# move these files into the project root is  probably best
source("~/Projects/Emulator/src/libRbind/EmuRbind.R")
source("~/Projects/Emulator/src/libRbind/multivar.R")

# defaults to the particle specta
sample <- gen.samp()

# do the pca decomposition
pca.decomp <- emu.pca(sample, nr=3) 

MAKETHETAS=FALSE

if(MAKETHETAS){
  thetas.est <- estimateThetas(pca.decomp)
  # output the thetas
  write.table(thetas.est, "theta-table.dat", row.names=FALSE, col.names=FALSE)
} else {
  thetas.est <- read.table("theta-table.dat")
}

r <- reconCurveAtPoint(c(0.15, 0.87, 0.78, 3.2, 1.64, 0.48, 0.19), thetas.est, pca.decomp)

plot(r$tvalues, r$mean, type="l", col="blue")
conf <- 1.956*sqrt(r$var)
lines(r$tvalues, r$mean+conf, col="red")
lines(r$tvalues, r$mean-conf, col="red")
