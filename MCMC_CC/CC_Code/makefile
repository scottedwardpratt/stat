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
INC= -I${GSLPATH}/include -I include -I ${INSTALLDIR}/include 
LIBS= -L${GSLPATH}/lib -L${INSTALLDIR}/lib

MCMC_HFILES = mcmc.h\
parameter.h\
distribution.h\
emulator.h

MCMC_OFILES = distribution.o\
emulator.o\
likelihood.o\
mcmc.o\
parameter.o\
prior.o\
proposal.o


SUBDIRS= build/mcmc\
lib\
include\

all: mcmcdirs lib/libmcmc.a install run
	
install: lib/libmcmc.a
	cp -f lib/libmcmc.a ${INSTALLDIR}/lib/;
	cp -f ${addprefix include/, ${MCMC_HFILES}} ${INSTALLDIR}/include;

mcmcdirs:
	mkdir -p ${SUBDIRS}
	
lib/libmcmc.a: ${addprefix include/, ${MCMC_HFILES}} ${addprefix build/mcmc/, ${MCMC_OFILES}}
	ar -ru lib/libmcmc.a ${addprefix build/mcmc/, ${MCMC_OFILES}};
	

include/%.h: src/%.h
	cp -f $< $@

build/mcmc/%.o: src/%.cc
	${CPP} -c ${OPT} ${INC} $< -o $@
	
clean:
	rm -f mcmc;
	
	for dir in ${SUBDIRS}; do\
		rm -f $$dir/*.o;\
		rm -f $$dir/*.h;\
		rm -f $$dir/*.a;\
	done
		
uninstall:
	rm  -f ${INSTALLDIR}/lib/libmcmc.a;
	rm  -f ${addprefix ${INSTALLDIR}/include/, ${MCMC_HFILES}}

run: mcmcrun.cc
	${CPP} ${OPT} ${INC} -o mcmc $< -lcoralutils -lgsl -lmcmc -lgslcblas ${LIBS}