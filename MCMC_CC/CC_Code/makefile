#enter make, then make install, make clean, make uninstall
#other commands: make coralutils, make coral, make xgraph, make install_coralutils, make install_coral, make install_xgraph
include ../makefile_defs.mk
#the variables below were set in the include statement above
CPP=${MADAI_CPP}
#e.g. /usr/bin/g++
OPT=${MADAI_CFLAGS}
#e.g. -O3
INSTALLDIR=${MADAI_INSTALLDIR}
#e.g. /User/scottepratt/local
GSLPATH=${MADAI_GSLPATH}
#e.g. /usr/local
#################################################
INC=-I ${GSLPATH}/include -I include -I $INSTALLDIR/include

MCMC_HFILES = mcmc.h\
parameter.h\
distribution.h

MCMC_OFILES = build/mcmc/proposal.o\
build/mcmc/mcmc.o


all: mcmcdirs lib/mcmc.a
	
mcmcdirs:
	mkdir -p build/mcmc lib include/
	
lib/mcmc.a: h_files ${MCMC_OFILES}
	mkdir -p lib;
	ar -ru lib/libmcmc.a ${MCMC_OFILES};
	
h_files:
	cp ${MCMC_HFILES} include/

build/mcmc/proposal.o: h_files proposal.cc
	${CPP} -c ${OPT} ${INC} proposal.cc -o build/mcmc/proposal.o	
	
build/mcmc/mcmc.o: h_files mcmc.cc
	${CPP} -c ${OPT} ${INC} mcmc.cc -o build/mcmc/mcmc.o

	