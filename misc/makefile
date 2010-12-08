#!/usr/bin/bash

CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}

all: ${MADAI_INCLUDE}/qualifier.h ${MADAI_BIN}/bmaker

${MADAI_INCLUDE}/qualifier.h: qualifier.h
	cp -f qualifier.h ${MADAI_INCLUDE}/qualifier.h
	
${MADAI_BIN}/bmaker: bmaker.cc
	${CPP} ${OPT} -o ${MADAI_BIN}/bmaker bmaker.cc

clean:
	rm ${MADAI_INCLUDE}/qualifier.h ${MADAI_BIN}/bmaker
