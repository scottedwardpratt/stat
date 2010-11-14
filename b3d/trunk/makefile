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
src/b3d/action.cc\
src/b3d/action_perform.cc\
src/b3d/action_perform_activate.cc\
src/b3d/action_perform_collide.cc\
src/b3d/action_perform_decay.cc\
src/b3d/action_perform_exitcell.cc\
src/b3d/action_perform_vizwrite.cc\
src/b3d/findcollision.cc\
src/b3d/hydrotob3d.cc\
src/b3d/part.cc\
src/b3d/resonances.cc

export ANALYZE_HFILES = src/b3d/b3d.h\
src/analyze/analyze.h

export ANALYZE_CCFILES = src/analyze/analyze.cc\
src/analyze/anal_balance.cc\
src/analyze/anal_gamma.cc\
src/analyze/anal_spectra.cc\
src/analyze/anal_v2.cc

CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}

all: b3d analyze

clean:
	rm ${MADAI_LIB}/libb3d.a ${MADAI_LIB}/libb3danalyze.a

b3d: headers ${MADAI_LIB}/libb3d.a
	for filename in ${B3D_HFILES};\
	do cp $${filename} ${MADAI_INCLUDE}/ ; done;

analyze: headers ${MADAI_LIB}/libb3danalyze.a
	for filename in ${ANALYZE_HFILES};\
	do cp $${filename} ${MADAI_INCLUDE}/ ; done;

headers: ${B3D_HFILES}
	for filename in ${B3D_HFILES};\
	do cp $${filename} include/ ; done;

${MADAI_LIB}/libb3d.a : ${B3D_HFILES} ${B3D_CCFILES}
	rm *.o build/b3d/*.o;\
	${CPP} -c -I ${MADAI_GSLPATH}/include -I${MADAI_HOME}/rhic/CorAL/trunk/include -I include -I ${MADAI_HDF5_HOME}/include -O ${OPT}  ${B3D_CCFILES} ;\
	mv *.o build/b3d/;\
	ar -ru lib/libb3d.a build/b3d/*.o;\
	mv -f lib/libb3d.a ${MADAI_LIB}/
	
${MADAI_LIB}/libb3danalyze.a : ${ANALYZE_HFILES} ${ANALYZE_CCFILES}
	rm *.o build/analyze/*.o;\
	${CPP} -c -I ${MADAI_GSLPATH}/include -I${MADAI_INCLUDE} -I ${MADAI_HDF5_HOME}/include -O ${OPT} ${ANALYZE_CCFILES};\
	mv *.o build/analyze/;\
	ar -ru lib/libanalyze.a build/analyze/*.o;\
	mv -f lib/libanalyze.a ${MADAI_LIB}/


	



