OPT=${MADAI_CFLAGS}
<<<<<<< HEAD
=======
CH=${CORAL_HOME}
CORALINCLUDE=${MADAI_INCLUDE}
CORALLIB=${MADAI_LIB}
GSLINC=${GSL_HOME}/include
GSLLIB=${GSL_HOME}/lib
>>>>>>> 46f911ce44427f0f5d6c3487a1669fb2796c2277

#HPM = -L{GSLLIB} -I{GSLINC} -L{CORALLIB} -lcoralutils 
# libraries needed for parameter map
HPM = -I${MADAI_INCLUDE} -I${MADAI_GSLPATH}/include -L${MADAI_GSLPATH}/lib -L${MADAI_LIB} -lcoralutils

clean:
	-rm -f *.o

CEos.o : CEos.cpp CEos.h def.h
	${MADAI_CPP} ${OPT} -I./ ${HPM} -c CEos.cpp -o CEos.o

CCell.o : CCell.cpp CCell.h CEos.o def.h
	${MADAI_CPP} ${OPT} -I./ ${HPM} -c CCell.cpp -o CCell.o 

CMesh.o : CMesh.cpp CMesh.h CCell.o CEos.o
	${MADAI_CPP} ${OPT} -I./ ${HPM} -c CMesh.cpp -o CMesh.o

isHydro3 : CCell.o CMesh.o CEos.o def.h main.cpp
	${MADAI_CPP} main.cpp -o isHydro3 -I./ ${OPT} ${HPM} CCell.o CMesh.o CEos.o

eos_test : CEos.o def.h eos_test.cpp
	${MADAI_CPP} -o eos_test -I./ ${OPT} ${HPM} eos_test.cpp CEos.o
