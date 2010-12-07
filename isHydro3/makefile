OPT=${MADAI_CFLAGS}

CH=${CORAL_HOME}
CORALINCLUDE=${MADAI_INCLUDE}
CORALLIB=${MADAI_LIB}
GSLINC=${GSL_HOME}/include
GSLLIB=${GSL_HOME}/lib

export HYDRO_HEADERS = CHydro.h CEos.h CCell.h CMesh.h hydroDef.h

# libraries needed for parameter map
HPM = -I${MADAI_INCLUDE} -I${MADAI_GSLPATH}/include -L${MADAI_GSLPATH}/lib -L${MADAI_LIB}

isHydro3 : ${MADAI_LIB}/libishydro3.a src/hydroDef.h main.cpp
	${MADAI_CPP} main.cpp -o isHydro3 ${OPT} ${HPM} -lishydro3 -lcoralutils

clean:
	rm -f build/*.o;\
	rm isHydro3\

build/CEos.o : src/CEos.cpp src/CEos.h src/hydroDef.h
	${MADAI_CPP} ${OPT} ${HPM} -c src/CEos.cpp -o build/CEos.o

build/CCell.o : src/CCell.cpp src/CCell.h build/CEos.o src/hydroDef.h
	${MADAI_CPP} ${OPT} ${HPM} -c src/CCell.cpp -o build/CCell.o 

build/CMesh.o : src/CMesh.cpp src/CMesh.h build/CCell.o build/CEos.o
	${MADAI_CPP} ${OPT} ${HPM} -c src/CMesh.cpp -o build/CMesh.o

build/CHydro.o : src/CHydro.cpp src/CHydro.h build/CMesh.o build/CCell.o build/CEos.o
	${MADAI_CPP} ${OPT} ${HPM} -c src/CHydro.cpp -o build/CHydro.o

${MADAI_LIB}/libishydro3.a : build/CHydro.o build/CMesh.o build/CCell.o build/CEos.o src/hydroDef.h
	ar -ru ${MADAI_LIB}/libishydro3.a build/*.o;\
	for filename in ${HYDRO_HEADERS};\
	do cp src/$${filename} ${MADAI_INCLUDE}/; done;\
