LHS sampler
ccs aug-10, cec24@phy.duke.edu

summary: for a model with a given number of parameters (nparams) produces a given number of 
samples (nruns) from an improved LHS design. The samples are all generated on the unit hyper-cube and
printed to the stdout.

requirements: linux/unix(probably), R, a shell

setup: install R (if you don't have it, you can test this by trying to start an R session)
			 install the lhs library, the R code will try and do this for you This can be done manually
			 by typing:
			
			 install.packages(c("lhs"))
			 
			 into a running R session. A little window will open prompting you to select a mirror then the 
			 package will download and install and be ready to use.

usage: sampler.sh takes 2 arguments:
			 			  npars -> the number of variables or parameters or dimensions to sample
							nruns -> the number of repeat samples to generate

			 eg: ./sampler.sh 3 10 
			 will generate a sample 10 vectors of 3 parameters each. 
			 										
contents: genSamples.R -> the r code which hooks into the lhs lib and generates the samples
					readme.txt -> the readme, you should look at this first
					sampler.sh -> shell wrapper, this is what you should run.
