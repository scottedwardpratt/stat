## same stuff as emuLum but with an emphasis on
## comparing the matern and power-exp cov fns
##
## this would ideally turn into a full integration
## over all potential cov-fns and regression orders


source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings

arch <- system("uname -s", intern=TRUE)
if(is.loaded("callEstimate") == FALSE){
  libNAME <- "~/local/lib/libRbind/libRBIND"
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

source("~/local/include/libRbind/emuOverDesign.R") # functions for running the emulator over the design in high dims
source("~/local/include/libRbind/implausOverDesign.R") # functions for computing the implausibility
source("~/local/include/libRbind/testEst.R") # for testing the estimation of thetas

# the model path for the first wave
modelOutputFile <- "./wave-1/lum_fun_outputs.dat"
modelData <- as.matrix(read.table(modelOutputFile))
nbins <- dim(modelData)[2]
nruns <- dim(modelData)[1]

# the design file for wave 1 (need to label this better)
designFile <- "./design/design_sorted_wave_1.dat"
desNames <- c("Zr", "Fescp", "Fbary")
nparams <- length(desNames)
designDefaults <- c(10, 50, 0.05)
# have to correct the design data, it currently includes an "unrun" run # 46 (i think)
designData <- as.matrix(read.table(designFile, col.names=desNames))[1:nruns,]


expDataFile <- "./lum_fun_observed.dat"
expDataIn <- as.matrix(read.table(expDataFile))
expData <- list( obsValue=expDataIn[2,], obsError=expDataIn[3,])


rebuild <- 1

if(rebuild == 1){
                                        # do the estimation for the 3 different cov fns
  estimResult.power <- doCombEstimation(cov.fn=1)
  save(estimResult.power, file="thetas-comb-power.dat")
  
  estimResult.matern32 <- doCombEstimation(cov.fn=2)
  save(estimResult.matern32, file="thetas-comb-matern32.dat")
  
  estimResult.matern52 <- doCombEstimation(cov.fn=3)
  save(estimResult.matern52, file="thetas-comb-matern52.dat")
}
