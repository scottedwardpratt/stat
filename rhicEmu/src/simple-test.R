## simple-test
## ccs, cec24@phy.duke.edu 2011-March/April (MSU)
##
## removes some number nhold from the training data given, 
## estimates hyperparams on the reduced set and then compares results
##
## mainly useful to get some idea if your current setup is working at all.
## this doesn't let you test a specific set of hyper-params unless you have already
## set aside some representative points to test the model on.

source("gen-samp.R") # load the sample set
source("emu-pca.R") # load the pca fns
# note, need to change this to .so or whatever for deployment
dyn.load("~/local/lib/libRBIND.dylib") #
# move these files into the project root is  probably best
source("~/local/include/libRbind/EmuRbind.R")

#read nhold from env
# this file reads the given sample set and generates a pca.decomp for it and the thetas
cargs <- Sys.getenv(c('infile', 'errfile', 'paramfile', 'nhold'))
inSamples <- cargs[1]
inSamplesErrs <- cargs[2]
inParams <- cargs[3]
nhold <- as.numeric(cargs[4])

cat("#witholding ", nhold, "points\n")

des.full <- load.design(inParams)
sample.full <- load.sample(inSamples, inSamplesErrs)

# now we want to remove some nhold points from the sample and the design
#  sequential design points are not close to each other so its ok to do it this way
des.cut <- list(nvars=des.full$nvars, nreps=(des.full$nreps-nhold), des=des.full$des[-(1:nhold),])
# the nhold held back points
des.held <- list(nvars=des.full$nvars, nreps=nhold, des=des.full$des[1:nhold,])

sample.cut <- list(nreps=(sample.full$nreps-nhold), ntps=sample.full$ntps,
                   y=sample.full$y[-(1:nhold),], means=sample.full$means,
                   sds=sample.full$sds, errs=sample.full$errs[-(1:nhold),])
# the nhold held back points
sample.held <- list(nreps=(nhold), ntps=sample.full$ntps,
                   y=sample.full$y[(1:nhold),], means=sample.full$means,
                   sds=sample.full$sds, errs=sample.full$errs[(1:nhold),])

# now we do estimation on the withheld sample
tvec <- seq(1, sample$ntps)
samplePCA.cut <- list(t=tvec, y=t(sample.cut$y), des=t(des.cut$des))

# do the pca decomp
pca.decomp.cut <- emu.pca(samplePCA.cut, nr=3, silent=TRUE)
save(pca.decomp.cut, file="PcaCutTemp.dat")

# now estimate the thetas
MAKETHETAS=FALSE
if(MAKETHETAS){
  thetas.est <- estimateThetas(pca.decomp.cut)
  write.table(thetas.est, "theta-held.dat", row.names=FALSE, col.names=FALSE)
} else {
  thetas.est <- read.table("theta-held.dat")
}

# now we can compute the values at each point
rvalues <- vector("list", npoints)

errors <- matrix(NA, nrow=nhold, ncol=sample$ntps)
vars <- matrix(NA, nrow=nhold, ncol=sample$ntps)

for(i in 1:nhold){
  # use the design points we held back to generate test values
  rvalues[[i]] <- reconCurveAtPoint(des.held$des[i,], thetas.est, pca.decomp.cut)
  # undo the centering and scaling we did in gen.sample
  # with the stored data from sample
  rvalues[[i]]$mean <- rvalues[[i]]$mean * sample$sds + sample$means
  rvalues[[i]]$var <- rvalues[[i]]$var * sample$sds
  #errors[i,] <- (sample.held$y[i,] - rvalues[[i]]$mean)  / sqrt(rvalues[[i]]$var)
  errors[i,] <- abs(rvalues[[i]]$mean - (sample.held$y[i,]*sample.held$sds + sample.held$means) )/ sqrt(rvalues[[i]]$var)
  vars[i,] <- rvalues[[i]]$var
}

write.table(errors, "nhold-errors.dat", row.names=FALSE, col.names=FALSE)
write.table(vars, "nhold-vars.dat", row.names=FALSE, col.names=FALSE)



                   





