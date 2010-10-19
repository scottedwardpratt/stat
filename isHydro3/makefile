OPT=${MADAI_CFLAGS}

CH=${CORAL_HOME}
CORALINCLUDE=${MADAI_INCLUDE}
CORALLIB=${MADAI_LIB}
GSLINC=${GSL_HOME}/include
GSLLIB=${GSL_HOME}/lib

# libraries needed for parameter map
HPM = -I${MADAI_INCLUDE} -I${MADAI_GSLPATH}/include -L${MADAI_GSLPATH}/lib -L${MADAI_LIB}

isHydro3 : ${MADAI_LIB}libishydro3.a def.h main.cpp
	${MADAI_CPP} main.cpp -o isHydro3 ${OPT} ${HPM} -lishydro3 -lcoralutils

clean:
	-rm -f *.o

CEos.o : CEos.cpp CEos.h def.h
	${MADAI_CPP} ${OPT} ${HPM} -c CEos.cpp -o CEos.o

CCell.o : CCell.cpp CCell.h CEos.o def.h
	${MADAI_CPP} ${OPT} ${HPM} -c CCell.cpp -o CCell.o 

CMesh.o : CMesh.cpp CMesh.h CCell.o CEos.o
	${MADAI_CPP} ${OPT} ${HPM} -c CMesh.cpp -o CMesh.o

CHydro.o : CHydro.cpp CHydro.h CMesh.o CCell.o CEos.o
	${MADAI_CPP} ${OPT} ${HPM} -c CHydro.cpp -o CHydro.o

${MADAI_LIB}libishydro3.a : CHydro.o CMesh.o CCell.o CEos.o
	ar -ru libishydro3.a *.o;\
	rm *.o;\
	mv -f libishydro3.a ${MADAI_LIB}

eos_test : CEos.o def.h eos_test.cpp
	${MADAI_CPP} -o eos_test -I./ ${OPT} ${HPM} eos_test.cpp CEos.o
