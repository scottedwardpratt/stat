#!/usr/bin/bash

INC=-I${MADAI_GSLPATH}/include -I${MADAI_INClUDE}/

all : ${MADAI_LIB}/libstat.a

${MADAI_LIB}/libstat.a : build/pca.o build/qualifier.o
	ar -ru lib/libstat.a build/pca.o build/qualifier.o;\
	mv -f lib/libstat.a ${MADAI_LIB}/
	
build/pca.o : src/pca.cc ${MADAI_INCLUDE}/pca.h ${MADAI_INCLUDE}/qualifier.h
	${MADAI_CPP} -c src/pca.cc ${MADAI_CFLAGS} -I${MADAI_GSLPATH}/include -I${MADAI_INCLUDE} -o build/pca.o

${MADAI_INCLUDE}/pca.h : include/pca.h
	cp include/pca.h ${MADAI_INCLUDE}/

include/pca.h : src/pca.h
	cp src/pca.h include/pca.h

build/qualifier.o : src/qualifier.cc ${MADAI_INCLUDE}/qualifier.h ${MADAI_INCLUDE}/pca.h
	${MADAI_CPP} -c src/qualifier.cc ${MADAI_CFLAGS} -I${MADAI_GSLPATH}/include -I${MADAI_INCLUDE} -o build/qualifier.o

${MADAI_INCLUDE}/qualifier.h : include/qualifier.h
	cp include/qualifier.h ${MADAI_INCLUDE}/

include/qualifier.h : src/qualifier.h
	cp src/qualifier.h include/qualifier.h




