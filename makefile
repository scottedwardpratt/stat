# enter make, make parmaker, make install, make clean, make uninstall
# local executables can be made with make runHydro or make multiEos
include ../makefile_defs.mk
#########################################################################
#HOPEFULLY, everything below need not be touched

INC=-I include -I${MADAI_GSLPATH}/include -I${MADAI_INSTALLDIR}/include -I${MADAI_HDF5_HOME}/include

all : lib/libstat.a parmaker

statdirs :
	mkdir -p lib;\
	mkdir -p include;\
	mkdir -p bin;\
	mkdir -p build;\
	mkdir -p ${MADAI_INSTALLDIR}/lib;\
	mkdir -p ${MADAI_INSTALLDIR}/include;\
	mkdir -p ${MADAI_INSTALLDIR}/bin

lib/libstat.a : build/pca.o build/qualifier.o
	rm lib/libstat.a;\
	ar -ru lib/libstat.a build/pca.o build/qualifier.o

build/pca.o : pca-src/pca.cc include/pca.h include/qualifier.h
	${MADAI_CPP} -c pca-src/pca.cc ${MADAI_CFLAGS} ${INC} -o build/pca.o

include/pca.h : pca-src/pca.h
	cp pca-src/pca.h include/pca.h

build/qualifier.o : pca-src/qualifier.cc include/qualifier.h include/pca.h
	${MADAI_CPP} -c pca-src/qualifier.cc ${MADAI_CFLAGS} ${INC} -o build/qualifier.o

include/qualifier.h : pca-src/qualifier.h
	cp pca-src/qualifier.h include/qualifier.h

parmaker : statdirs bin/parmaker

bin/parmaker : latinhyper3/parmaker.cc
	${MADAI_CPP} ${MADAI_CFLAGS} latinhyper3/parmaker.cc -o bin/parmaker

#########################

install : ${MADAI_INSTALLDIR}/lib/libstat.a ${MADAI_INSTALLDIR}/include/pca.h ${MADAI_INSTALLDIR}/include/qualifier.h ${MADAI_INSTALLDIR}/bin/parmaker ${MADAI_INSTALLDIR}/progdata/genSamples.R ${MADAI_INSTALLDIR}/bin/latin3.sh ${MADAI_INSTALLDIR}/bin/latin3.sh

${MADAI_INSTALLDIR}/lib/libstat.a : lib/libstat.a
	rm -f ${MADAI_INSTALLDIR}/lib/libstat.a;
	cp -f lib/libstat.a ${MADAI_INSTALLDIR}/lib/

${MADAI_INSTALLDIR}/include/pca.h : include/pca.h
	cp include/pca.h ${MADAI_INSTALLDIR}/include

${MADAI_INSTALLDIR}/include/qualifier.h : include/qualifier.h
	cp include/qualifier.h ${MADAI_INSTALLDIR}/include

${MADAI_INSTALLDIR}/bin/parmaker : bin/parmaker
	cp -f bin/parmaker ${MADAI_INSTALLDIR}/bin/

${MADAI_INSTALLDIR}/progdata/genSamples.R : latinhyper3/genSamples.R
	cp -f latinhyper3/genSamples.R ${MADAI_INSTALLDIR}/progdata/

${MADAI_INSTALLDIR}/bin/latin3.sh : latinhyper3/latin3.sh
	cp -f latinhyper3/latin3.sh ${MADAI_INSTALLDIR}/bin/

clean :
	rm -f include/* build/* lib/* bin/*

uninstall :
	rm -f ${MADAI_INSTALLDIR}/lib/libstat.a ${MADAI_INSTALLDIR}/include/pca.h ${MADAI_INSTALLDIR}/include/qualifier.h ${MADAI_INSTALLDIR}/bin/parmaker ${MADAI_INSTALLDIR}/progdata/genSamples.R ${MADAI_INSTALLDIR}/bin/latin3.sh ${MADAI_INSTALLDIR}/bin/latin3.sh
