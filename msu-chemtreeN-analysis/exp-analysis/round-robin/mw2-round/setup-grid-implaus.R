## particular variables for the mw2 analysis
## we're going to compare against mw2 to start with
mwString = "mw2"
mwString.Up = toupper(mwString)

compareStrings <- c("mw1", "mw2", "mw3", "mw6")
ncompare <- length(compareStrings)

## do this for everything...
for(index in 1:ncompare){
  mwStringCompare = compareStrings[index]
  mwStringCompare.Up = toupper(mwStringCompare)

  cat("# mwCompare: " , mwStringCompare, "\n")
  source("../gridCompareImplaus.R")
}
