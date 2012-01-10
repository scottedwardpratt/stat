## the path to the design file
designFile <- "./design/stats.param-combined-s"
## change this to pick which cclass we'll look at
cclass <- "0to10"
## path to the model data
dataPath <- "./model-data"

## base name of the analysis (could be changed automatically)
analysisBaseName <- "spectra-combined-"
## modelData
modelDataFile <- paste(dataPath, "/", analysisBaseName , cclass, ".dat", sep="")
## the experimental data
defaultFile <- paste(dataPath, "/", analysisBaseName, cclass, "-default.dat", sep="")
## the per run errors (if we care about these)
errorFile <- paste(dataPath, "/", analysisBaseName, cclass, "-errors.dat", sep="")


designFieldNames <- c("GLAUBER_WNBIN_RATIO", "GLAUBER_K_TAU",
                      "HYDRO_INIT_NS", "HYDRO_INIT_FLOW", "HYDRO_SVRATIO",
                      "EQOFST_BUMP_HEIGHT")
