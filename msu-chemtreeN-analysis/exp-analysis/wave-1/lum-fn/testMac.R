## ccs, there seems to be some issue with doing a source on this file when rebuild
## is set to 1, on my os-x 10.7 mbp at least. Is this an R issue?

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

# do we want to rebuild the hyperparameters or not? 1 for yes
buffer <- paste("./thetas-comb-lum.dat")

## lets try another external command instead of doCombEstimation
## estimResult1 <- doCombEstimation(fixNugget=NULL)
## cat("did run 1\n")
## estimResult2 <- doCombEstimation(fixNugget=NULL)
## cat("did run 2\n")

# plot the implausibility against dims 1 and 2
impStepPlots <- function(){
  for(i in 1:nbins){
    buffer <- paste("./images/implaus-stepped-cent-bin-",i,".pdf", sep="")
    cat("implausStepping: ", buffer, "\n")
    pdf(buffer)
    stepPlotDimensionImplaus(i,1,3,2, estim.result=estimResult, exp.data=expData)
    dev.off()
  }
}


load(buffer)
impStepPlots()
cat("did plots once\n")
impStepPlots()
# and we crash, wtf R
cat("did plots again\n")



