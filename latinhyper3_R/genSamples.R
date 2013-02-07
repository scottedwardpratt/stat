## LHS sampler
## ccs aug-10, cec24@phy.duke.edu
## RELEASED UNDER NO WARRANTY

retval <- require("lhs", quietly=TRUE)
options(echo=FALSE)

## test for the library
if(retval == FALSE){
  print("installing the lhs library")
  install.packages(c("lhs"))
}

generateSamples <- function(npars, nruns){
  ## we'll use the improvedLHS fn from the LHS library
  ## there are a few other generators you could use, but since
  ## a nearly max-min distribution is pretty ideal this is a good choice.
  data <- improvedLHS(nruns, npars)
  write.table(data, row.names=FALSE, col.names=FALSE, sep="\t")
}

cargs <- Sys.getenv(c('npars', 'nruns'))
nparameters <- as.numeric(cargs[1])
nruns <- as.numeric(cargs[2])
generateSamples(nparameters, nruns)
