OPT=${MADAI_CFLAGS}
HPM = -I${MADAI_INCLUDE} -I${MADAI_GSLPATH}/include -L${MADAI_GSLPATH}/lib -L${MADAI_LIB}
thisdir = ${PWD}

all: b3d analyze isHydro3

isHydro3 :  hydromain.cc hydroLibs
	${MADAI_CPP} -o isHydro3 hydromain.cc ${OPT} ${HPM} -lcoralutils -lishydro3

libs :
	make -C ${MADAI_HOME}/rhic/b3d/trunk -f ${MADAI_HOME}/rhic/b3d/trunk/makefile all

hydroLibs :
	make -C ${MADAI_HOME}/rhic/isHydro3 -f ${MADAI_HOME}/rhic/isHydro3/makefile isHydro3

eoscalc : eoscalc.cc libs
	${MADAI_CPP} -o eoscalc eoscalc.cc  ${OPT} ${HPM} -lcoralutils -lb3d -lgsl

b3d : b3dmain.cc libs
	${MADAI_CPP} -o b3d b3dmain.cc ${OPT} ${HPM} -lb3d -lcoralutils -lgsl -lgslcblas -lhdf5 -lhdf5_cpp -lhdf5_hl -lhdf5_hl_cpp -lsz -lz
	
analyze : analyzemain.cc libs
		${MADAI_CPP} -o analyze analyzemain.cc ${OPT} ${HPM} -lcoralutils -lanalyze -lgsl -lgslcblas -lhdf5 -lhdf5_cpp -lhdf5_hl -lhdf5_hl_cpp -lsz -lz