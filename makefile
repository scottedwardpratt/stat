# enter make, make parmaker, make install, make clean, make uninstall
# local executables can be made with make runHydro or make multiEos
include ../rhic/makefile_defs.mk
### the variables below are set by the include statement above
CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}
CORAL_INCLUDE=${INSTALLDIR}/include
#e.g. /Users/scottepratt/git -- 
CORAL_HOME=${MADAI_HOME}/rhic/CorAL/trunk
#e.g. ${MADAI_HOME}/rhic/CorAL/trunk
INSTALLDIR=${MADAI_INSTALLDIR}
#e.g. /Users/scottepratt/local
GSLPATH=${MADAI_GSLPATH}
#e.g. /opt/local
#########################################################################
#HOPEFULLY, everything below need not be touched

INC=-I include -I${GSLPATH}/include -I${INSTALLDIR}/include

all : lib/libstat.a parmaker

statdirs :
	mkdir -p lib;\
	mkdir -p include;\
	mkdir -p bin;\
	mkdir -p build;\
	mkdir -p ${INSTALLDIR}/lib;\
	mkdir -p ${INSTALLDIR}/include;\
	mkdir -p ${INSTALLDIR}/bin

lib/libstat.a : statdirs build/pca.o build/qualifier.o
	ar -ru lib/libstat.a build/pca.o build/qualifier.o

build/pca.o : src/pca.cc include/pca.h include/qualifier.h
	${CPP} -c src/pca.cc ${OPT} ${INC} -o build/pca.o

include/pca.h : src/pca.h
	cp src/pca.h include/pca.h

build/qualifier.o : src/qualifier.cc include/qualifier.h include/pca.h
	${CPP} -c src/qualifier.cc ${OPT} -I${INC} -o build/qualifier.o

include/qualifier.h : src/qualifier.h
	cp src/qualifier.h include/qualifier.h

parmaker : statdirs bin/parmaker

bin/parmaker : latinhyper3/parmaker.cc
	${CPP} ${OPT} latinhyper3/parmaker.cc -o bin/parmaker


#########################

install : ${INSTALLDIR}/lib/libstat.a ${INSTALLDIR}/include/pca.h ${INSTALLDIR}/include/qualifier.h ${INSTALLDIR}/bin/parmaker ${INSTALLDIR}/progdata/genSamples.R ${INSTALLDIR}/bin/latin3.sh ${INSTALLDIR}/bin/latin3.sh

${INSTALLDIR}/lib/libstat.a : lib/libstat.a
	cp -f lib/libstat.a ${INSTALLDIR}/lib/

${INSTALLDIR}/include/pca.h : include/pca.h
	cp include/pca.h ${INSTALLDIR}/include

${INSTALLDIR}/include/qualifier.h : include/qualifier.h
	cp include/qualifier.h ${INSTALLDIR}/include

${INSTALLDIR}/bin/parmaker : bin/parmaker
	cp -f bin/parmaker ${INSTALLDIR}/bin/

${INSTALLDIR}/progdata/genSamples.R : latinhyper3/genSamples.R
	cp -f latinhyper3/genSamples.R ${INSTALLDIR}/progdata/

${INSTALLDIR}/bin/latin3.sh : latinhyper3/latin3.sh
	cp -f latinhyper3/latin3.sh ${INSTALLDIR}/bin/

clean :
	rm -f include/* build/* lib/* bin/*

uninstall :
	rm -f ${INSTALLDIR}/lib/libstat.a ${INSTALLDIR}/include/pca.h ${INSTALLDIR}/include/qualifier.h ${INSTALLDIR}/bin/parmaker ${INSTALLDIR}/progdata/genSamples.R ${INSTALLDIR}/bin/latin3.sh ${INSTALLDIR}/bin/latin3.sh
