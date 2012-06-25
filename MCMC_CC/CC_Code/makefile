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
INC= -I${GSLPATH}/include -I include -I ${INSTALLDIR}/include -I ${INSTALLDIR}/include/EmuPlusPlus
LIBS= -L${GSLPATH}/lib -L${INSTALLDIR}/lib

MCMC_HFILES = mcmc.h\
parameterset.h\
distribution.h\
visualization.h\
quad.h
#emulator.h\ #This is no longer used because we use the interactive emulator

MCMC_OFILES = distribution.o\
quad.o\
likelihood_rhic.o\
mcmc.o\
parameterset.o\
prior_rhic.o\
proposal.o\
visualization.o\
likelihood_rhic_pca.o\
likelihood_cosmo.o\
likelihood_test.o\
prior_rhic_pca.o\
prior_cosmo.o\
prior_test.o\
#emulator.o\ #This is no longer used because we use the interactive emulator


SUBDIRS= build\
lib\
include\

all: mcmcdirs lib/libmcmc.a

temp: all install run
	
install: lib/libmcmc.a
	cp -f lib/libmcmc.a ${INSTALLDIR}/lib/;
	cp -f ${addprefix include/, ${MCMC_HFILES}} ${INSTALLDIR}/include;

mcmcdirs:
	mkdir -p ${SUBDIRS}
	
lib/libmcmc.a: ${addprefix include/, ${MCMC_HFILES}} ${addprefix build/, ${MCMC_OFILES}}
	rm -f lib/libmcmc.a;
	ar -ru lib/libmcmc.a ${addprefix build/, ${MCMC_OFILES}};


include/%.h: src/Headers/%.h
	cp -f $< $@

build/likelihood_rhic.o: src/RHIC/likelihood_rhic.cc
	${CPP} -c ${OPT} ${INC} $< -o $@

build/likelihood_rhic_pca.o: src/RHIC_PCA/likelihood_rhic_pca.cc
	${CPP} -c ${OPT} ${INC} $< -o $@

build/likelihood_cosmo.o: src/Cosmosurvey/likelihood_cosmo.cc
	${CPP} -c ${OPT} ${INC} $< -o $@
	
build/likelihood_test.o: src/TEST/likelihood_test.cc
	${CPP} -c ${OPT} ${INC} $< -o $@
	
build/prior_rhic.o: src/RHIC/prior_rhic.cc
	${CPP} -c ${OPT} ${INC} $< -o $@

build/prior_rhic_pca.o: src/RHIC_PCA/prior_rhic_pca.cc
	${CPP} -c ${OPT} ${INC} $< -o $@

build/prior_cosmo.o: src/Cosmosurvey/prior_cosmo.cc
	${CPP} -c ${OPT} ${INC} $< -o $@
	
build/prior_test.o: src/TEST/prior_test.cc
	${CPP} -c ${OPT} ${INC} $< -o $@
	
build/%.o: src/Common_Files/%.cc
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
	${CPP} ${OPT} ${INC} -o mcmc $< -lcoralutils -lgsl -lmcmc -lgslcblas -lEmuPlusPlus ${LIBS}
	cp mcmc ~/local/include/;

test: emulatortest.cc
	${CPP} ${OPT} ${INC} -o emutest $< -lcoralutils -lgsl -lmcmc -lgslcblas ${LIBS}
