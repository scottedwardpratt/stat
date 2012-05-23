# enter make, make parmaker, make install, make clean, make uninstall
# local executables can be made with make runHydro or make multiEos
include ../../makefile_defs.mk
#########################################################################
#HOPEFULLY, everything below need not be touched

INC=-I include -I${MADAI_GSLPATH}/include -I${MADAI_INSTALLDIR}/include -I${MADAI_HDF5_HOME}/include -I${MADAI_GSLPATH}/include

LB=-L${MADAI_INSTALLDIR}/lib -L${MADAI_GSLPATH}/lib -lcoralutils -lgsl

all : lib/librhicstat.a parmaker_random parmaker_latin

statdirs :
	mkdir -p lib;\
	mkdir -p include;\
	mkdir -p bin;\
	mkdir -p build;\
	mkdir -p ${MADAI_INSTALLDIR}/lib;\
	mkdir -p ${MADAI_INSTALLDIR}/include;\
	mkdir -p ${MADAI_INSTALLDIR}/bin

lib/librhicstat.a : build/pca.o build/qualifier.o build/rhicstat.o build/calcsens.o build/quadfit.o
	rm -f lib/librhicstat.a;\
	ar -ru lib/librhicstat.a build/pca.o build/qualifier.o build/rhicstat.o build/calcsens.o build/quadfit.o

build/pca.o : src/pca.cc include/rhicstat.h include/qualifier.h
	${MADAI_CPP} -c src/pca.cc ${MADAI_CFLAGS} ${INC} -o build/pca.o
	
build/rhicstat.o : src/rhicstat.cc include/rhicstat.h include/qualifier.h
	${MADAI_CPP} -c src/rhicstat.cc ${MADAI_CFLAGS} ${INC} -o build/rhicstat.o

build/calcsens.o : src/calcsens.cc include/rhicstat.h include/qualifier.h
	${MADAI_CPP} -c src/calcsens.cc ${MADAI_CFLAGS} ${INC} -o build/calcsens.o

build/quadfit.o : src/quadfit.cc include/rhicstat.h include/qualifier.h
	${MADAI_CPP} -c src/quadfit.cc ${MADAI_CFLAGS} ${INC} -o build/quadfit.o

include/rhicstat.h : src/rhicstat.h
	cp src/rhicstat.h include/rhicstat.h

build/qualifier.o : src/qualifier.cc include/qualifier.h include/rhicstat.h
	${MADAI_CPP} -c src/qualifier.cc ${MADAI_CFLAGS} ${INC} -o build/qualifier.o

include/qualifier.h : src/qualifier.h
	cp src/qualifier.h include/qualifier.h

parmaker : statdirs bin/parmaker_latin bin/parmaker_random

bin/parmaker_latin : latinhyper3/parmaker_latin.cc
	${MADAI_CPP} ${MADAI_CFLAGS} latinhyper3/parmaker_latin.cc -o bin/parmaker_latin

bin/parmaker_random : src/parmaker_random.cc
	${MADAI_CPP} ${MADAI_CFLAGS} ${INC} ${LB} src/parmaker_random.cc -o bin/parmaker_random

#########################

install : ${MADAI_INSTALLDIR}/lib/librhicstat.a ${MADAI_INSTALLDIR}/include/rhicstat.h ${MADAI_INSTALLDIR}/include/qualifier.h ${MADAI_INSTALLDIR}/bin/parmaker_latin ${MADAI_INSTALLDIR}/bin/parmaker_random ${MADAI_INSTALLDIR}/progdata/genSamples.R ${MADAI_INSTALLDIR}/bin/latin3.sh ${MADAI_INSTALLDIR}/bin/latin3.sh

${MADAI_INSTALLDIR}/lib/librhicstat.a : lib/librhicstat.a
	rm -f ${MADAI_INSTALLDIR}/lib/librhicstat.a;
	cp -f lib/librhicstat.a ${MADAI_INSTALLDIR}/lib/


${MADAI_INSTALLDIR}/include/rhicstat.h : include/rhicstat.h
	cp include/rhicstat.h ${MADAI_INSTALLDIR}/include

${MADAI_INSTALLDIR}/include/qualifier.h : include/qualifier.h
	cp include/qualifier.h ${MADAI_INSTALLDIR}/include

${MADAI_INSTALLDIR}/bin/parmaker_random : bin/parmaker_random
	cp -f bin/parmaker_random ${MADAI_INSTALLDIR}/bin/

${MADAI_INSTALLDIR}/bin/parmaker_latin : bin/parmaker_latin
	cp -f bin/parmaker_latin ${MADAI_INSTALLDIR}/bin/


${MADAI_INSTALLDIR}/progdata/genSamples.R : latinhyper3/genSamples.R
	cp -f latinhyper3/genSamples.R ${MADAI_INSTALLDIR}/progdata/

${MADAI_INSTALLDIR}/bin/latin3.sh : latinhyper3/latin3.sh
	cp -f latinhyper3/latin3.sh ${MADAI_INSTALLDIR}/bin/

clean :
	rm -f include/* build/* lib/* bin/*

uninstall :
	rm -f ${MADAI_INSTALLDIR}/lib/librhicstat.a ${MADAI_INSTALLDIR}/include/rhicstat.h ${MADAI_INSTALLDIR}/include/qualifier.h ${MADAI_INSTALLDIR}/bin/parmaker_latin ${MADAI_INSTALLDIR}/bin/parmaker_random  ${MADAI_INSTALLDIR}/progdata/genSamples.R ${MADAI_INSTALLDIR}/bin/latin3.sh ${MADAI_INSTALLDIR}/bin/latin3.sh
