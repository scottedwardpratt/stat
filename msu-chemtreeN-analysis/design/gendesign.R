#######################################################
# gendesign  summer - 2011
# ccs cec24@phy.duke.edu
# 
# routines for generating a design for the madai-emulator
# applied to facundo's galaxy model called (?)
#
# see info.txt for more details, we currently have 3 params
# Zr \in [5,20] 10
# Fescp \in [0,100] 50
# Fbary \in [0,0.2] 0.05
#
# for fun, the final run should be the default values
######################################################

# load the Latin Hypercube Sampling pkg
# if you don't have it, install that bastard
retval <- require("lhs", quietly=TRUE)                    
#options(echo=FALSE)
if(retval == FALSE){
  print("#installing the lhs library")
  install.packages(c("lhs"))
}


nDesignPoints <- 45
nParams <- 3


paramRanges <- matrix(0, nrow=nParams, ncol=3)
# we ought to load this from the info file... (but i'm lazy)
paramRanges[1,] <- c(5,20,10)
paramRanges[2,] <- c(0,100,50)
paramRanges[3,] <- c(0,0.2,0.05)

# generate the design on the unit (hyper)cube 
design.unscaled <- maximinLHS(nDesignPoints, nParams)
# this is the matrix we'll scale up to fit on the parameter range
design.scaled <- matrix(0, nrow=nDesignPoints, ncol=nParams)

# scale the data
for(i in 1:nParams)
  design.scaled[,i] <- (design.unscaled[,i]*(paramRanges[i,2] - paramRanges[i,1]) + paramRanges[i,1])

# now add the default point
design.final <- rbind(design.scaled, paramRanges[,3])
filename <- "design.dat"
cat("# design file for galaxy-experiment\n", file=filename)
cat("# nDesignPoints ", nDesignPoints, " nParams: ", nParams , "\n", file=filename, append=TRUE)
cat("# default run is the final run listed", file=filename, append=TRUE)
cat("# min, max, default\n", file=filename, append=TRUE)

for(i in 1:nParams)
  cat("# ", paramRanges[i,1], " ", paramRanges[i,2], " ", paramRanges[i,3], "\n", file=filename, append=TRUE)

cat("# Zr\tFescp\tFbary\n", file=filename, append=TRUE)
write.table(design.final, file="design.dat", col.names=FALSE, row.names=FALSE, append=TRUE)



