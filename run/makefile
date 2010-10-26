OPT=${MADAI_CFLAGS}
HPM = -I${MADAI_INCLUDE} -I${MADAI_GSLPATH}/include -L${MADAI_GSLPATH}/lib -L${MADAI_LIB}
 
isHydro3 :  hydromain.cc
	${MADAI_CPP} -o isHydro3 hydromain.cc -I./ ${OPT} ${HPM} -lcoralutils -lishydro3

eoscalc : eoscalc.cc
	${MADAI_CPP} -o eoscalc eoscalc.cc  ${OPT} ${HPM} -lcoralutils -lb3d -lgsl