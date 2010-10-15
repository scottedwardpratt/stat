OPT=${CFLAGS}
CH=${CORAL_HOME}
CORALINCLUDE=${CH}/include
CORALLIB=${CH}/lib
GSLINC=${GSL_HOME}/include
GSLLIB=${GSL_HOME}/lib

#HPM = -L{GSLLIB} -I{GSLINC} -L{CORALLIB} -lcoralutils 
# libraries needed for parameter map
HPM = -I${CORALINCLUDE} -L${GSLLIB} -I${GSLINC} -L${CORALLIB} 

clean:
	-rm -f *.o

CEos.o : CEos.cpp CEos.h def.h
	g++ ${OPT} -I./ ${HPM} -c CEos.cpp -o CEos.o

CCell.o : CCell.cpp CCell.h CEos.o def.h
	g++ ${OPT} -I./ ${HPM} -c CCell.cpp -o CCell.o 

CMesh.o : CMesh.cpp CMesh.h CCell.o CEos.o
	g++ ${OPT} -I./ ${HPM} -c CMesh.cpp -o CMesh.o

isHydro3 : CCell.o CMesh.o CEos.o def.h main.cpp
	g++ -o isHydro3 -I./ ${OPT} ${HPM} -lcoralutils main.cpp CCell.o CMesh.o CEos.o

eos_test : CEos.o def.h eos_test.cpp
	g++ -o eos_test -I./ ${OPT} ${HPM} -lcoralutils eos_test.cpp CEos.o
