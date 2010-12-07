OPT=${MADAI_CFLAGS}

CH=${CORAL_HOME}
CORALINCLUDE=${MADAI_INCLUDE}
CORALLIB=${MADAI_LIB}
GSLINC=${GSL_HOME}/include
GSLLIB=${GSL_HOME}/lib

export HYDRO_HEADERS = src/CHydro.h src/CEos.h src/CCell.h src/CMesh.h src/hydroDef.h

# libraries needed for parameter map
HPM = -I${MADAI_INCLUDE} -I${MADAI_GSLPATH}/include -L${MADAI_GSLPATH}/lib -L${MADAI_LIB}

isHydro3 : ${MADAI_LIB}/libishydro3.a src/hydroDef.h main.cpp
	${MADAI_CPP} main.cpp -o isHydro3 ${OPT} ${HPM} -lishydro3 -lcoralutils

clean:
	rm -f build/*.o
	rm isHydro3

build/CEos.o : src/CEos.cpp src/CEos.h src/hydroDef.h
	${MADAI_CPP} ${OPT} ${HPM} -c src/CEos.cpp -o build/CEos.o

build/CCell.o : src/CCell.cpp src/CCell.h build/CEos.o
	${MADAI_CPP} ${OPT} ${HPM} -c src/CCell.cpp -o build/CCell.o 

build/CMesh.o : src/CMesh.cpp src/CMesh.h build/CCell.o build/CEos.o
	${MADAI_CPP} ${OPT} ${HPM} -c src/CMesh.cpp -o build/CMesh.o

build/CHydro.o : src/CHydro.cpp src/CHydro.h build/CMesh.o build/CCell.o build/CEos.o
	${MADAI_CPP} ${OPT} ${HPM} -c src/CHydro.cpp -o build/CHydro.o

${MADAI_LIB}/libishydro3.a : build/CHydro.o build/CMesh.o build/CCell.o build/CEos.o 
	ar -ru libishydro3.a build/*.o;\
	for filename in ${HYDRO_HEADERS};\
	do cp $${filename} ${MADAI_INCLUDE}/; done;\
	mv -f libishydro3.a ${MADAI_LIB}/
