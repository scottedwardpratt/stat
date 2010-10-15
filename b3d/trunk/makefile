#!/usr/bin/bash
export B3D_HFILES = src/b3d/b3d.h\
src/b3d/bjmaker.h\
src/b3d/resonances.h\
src/analyze/analyze.h

export B3D_CCFILES = src/b3d/b3d.cc\
src/b3d/bjmaker.cc\
src/b3d/cell.cc\
src/b3d/collide.cc\
src/b3d/decay.cc\
src/b3d/event.cc\
src/b3d/event_perform.cc\
src/b3d/event_perform_activate.cc\
src/b3d/event_perform_collide.cc\
src/b3d/event_perform_decay.cc\
src/b3d/event_perform_exitcell.cc\
src/b3d/findcollision.cc\
src/b3d/hydrotob3d.cc\
src/b3d/part.cc\
src/b3d/resonances.cc

CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}

b3d: headers ${MADAI_LIB}/libb3d.a ${MADAI_LIB}/libb3danalyze.a
	for filename in ${B3D_HFILES};\
		do cp $${filename} ${MADAI_INCLUDE}/ ; done;

headers: ${B3D_HFILES}
	for filename in ${B3D_HFILES};\
	do cp $${filename} include/ ; done;

${MADAI_LIB}/libb3d.a : ${B3D_HFILES} ${B3D_CCFILES}
	rm *.o build/b3d/*.o;\
	${CPP} -c -I ${MADAI_GSLPATH}/include -I${MADAI_HOME}/rhic/CorAL/trunk/include -I include -I ${MADAI_HDF5_HOME}/include -O ${OPT}  ${B3D_CCFILES};\
	mv *.o build/b3d/;\
	ar -ru lib/libb3d.a build/b3d/*.o;\
	mv -f lib/libb3d.a ${MADAI_LIB}/
	
${MADAI_LIB}/libb3danalyze.a : src/analyze/analyze.cc
	rm *.o build/analyze/*.o;\
	${CPP} -c -I ${MADAI_GSLPATH}/include -I${MADAI_INCLUDE} -I ${MADAI_HDF5_HOME}/include -O ${OPT} src/analyze/analyze.cc;\
	mv *.o build/analyze/;\
	ar -ru lib/libanalyze.a build/analyze/*.o;\
	mv -f lib/libanalyze.a ${MADAI_LIB}/


	



