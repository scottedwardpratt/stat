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

export B3D_OBJFILES = build/b3d/b3d.o\
build/b3d/bjmaker.o\
build/b3d/cell.o\
build/b3d/collide.o\
build/b3d/decay.o\
build/b3d/action.o\
build/b3d/action_perform.o\
build/b3d/action_perform_activate.o\
build/b3d/action_perform_collide.o\
build/b3d/action_perform_decay.o\
build/b3d/action_perform_exitcell.o\
build/b3d/action_perform_vizwrite.o\
build/b3d/findcollision.o\
build/b3d/hydrotob3d.o\
build/b3d/part.o\
build/b3d/resonances.o

export ANALYZE_HFILES = src/b3d/b3d.h\
src/analyze/analyze.h

export ANALYZE_CCFILES = src/analyze/analyze.cc\
src/analyze/anal_balance.cc\
src/analyze/anal_gamma.cc\
src/analyze/anal_spectra.cc\
src/analyze/anal_v2.cc

export ANALYZE_OBJFILES = build/analyze/analyze.o\
build/analyze/anal_balance.o\
build/analyze/anal_gamma.o\
build/analyze/anal_spectra.o\
build/analyze/anal_v2.o

CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}
INC=-I ${MADAI_GSLPATH}/include -I${MADAI_HOME}/rhic/CorAL/trunk/include -I include -I ${MADAI_HDF5_HOME}/include
LIB=-L${MADAI_GSLPATH}/lib -L${MADAI_LIB} -lcoralutils -lanalyze -lgsl -lhdf5 -lhdf5_cpp -lhdf5_hl -lhdf5_hl_cpp -lsz -lz

all: b3d analyze

clean:
	rm ${MADAI_LIB}/libb3d.a ${MADAI_LIB}/libb3danalyze.a build/b3d/*.o build/analyze/*.o

b3d: headers ${MADAI_LIB}/libb3d.a

analyze: headers ${MADAI_LIB}/libb3danalyze.a

headers: ${B3D_HFILES}
	for filename in ${B3D_HFILES};\
	do cp $${filename} ${MADAI_INCLUDE}/ ; done;

${MADAI_LIB}/libb3d.a : ${B3D_HFILES} ${B3D_OBJFILES}
	ar -ru lib/libb3d.a ${B3D_OBJFILES};\
	mv -f lib/libb3d.a ${MADAI_LIB}/

build/b3d/b3d.o : src/b3d/b3d.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/b3d.o src/b3d/b3d.cc
	
build/b3d/bjmaker.o : src/b3d/bjmaker.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/bjmaker.o src/b3d/bjmaker.cc

build/b3d/cell.o : src/b3d/cell.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/cell.o src/b3d/cell.cc

build/b3d/collide.o : src/b3d/collide.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/collide.o src/b3d/collide.cc

build/b3d/decay.o : src/b3d/decay.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/decay.o src/b3d/decay.cc

build/b3d/action.o : src/b3d/action.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/action.o src/b3d/action.cc

build/b3d/action_perform.o : src/b3d/action_perform.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/action_perform.o src/b3d/action_perform.cc

build/b3d/action_perform_activate.o : src/b3d/action_perform_activate.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/action_perform_activate.o src/b3d/action_perform_activate.cc

build/b3d/action_perform_collide.o : src/b3d/action_perform_collide.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/action_perform_collide.o src/b3d/action_perform_collide.cc

build/b3d/action_perform_decay.o : src/b3d/action_perform_decay.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/action_perform_decay.o src/b3d/action_perform_decay.cc

build/b3d/action_perform_exitcell.o : src/b3d/action_perform_exitcell.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/action_perform_exitcell.o src/b3d/action_perform_exitcell.cc
	
build/b3d/action_perform_vizwrite.o : src/b3d/action_perform_vizwrite.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/action_perform_vizwrite.o src/b3d/action_perform_vizwrite.cc
	
build/b3d/findcollision.o : src/b3d/findcollision.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/findcollision.o src/b3d/findcollision.cc
	
build/b3d/hydrotob3d.o : src/b3d/hydrotob3d.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/hydrotob3d.o src/b3d/hydrotob3d.cc
	
build/b3d/part.o : src/b3d/part.cc
	${CPP} -c ${OPT} ${INC} -o build/b3d/part.o src/b3d/part.cc
		
build/b3d/resonances.o : 	src/b3d/resonances.cc
		${CPP} -c ${OPT} ${INC} -o build/b3d/resonances.o src/b3d/resonances.cc

${MADAI_LIB}/libb3danalyze.a : ${ANALYZE_HFILES} ${ANALYZE_OBJFILES}
	ar -ru lib/libb3danalyze.a ${ANALYZE_OBJFILES};\
	mv -f lib/libb3danalyze.a ${MADAI_LIB}/

build/analyze/analyze.o : src/analyze/analyze.cc
	${CPP} -c ${OPT} ${INC} -o build/analyze/analyze.o src/analyze/analyze.cc
	
build/analyze/anal_balance.o : src/analyze/anal_balance.cc
	${CPP} -c ${OPT} ${INC} -o build/analyze/anal_balance.o src/analyze/anal_balance.cc
		
build/analyze/anal_gamma.o : src/analyze/anal_gamma.cc
	${CPP} -c ${OPT} ${INC} -o build/analyze/anal_gamma.o src/analyze/anal_gamma.cc

build/analyze/anal_spectra.o : src/analyze/anal_spectra.cc
	${CPP} -c ${OPT} ${INC} -o build/analyze/anal_spectra.o src/analyze/anal_spectra.cc

build/analyze/anal_v2.o : src/analyze/anal_v2.cc
		${CPP} -c ${OPT} ${INC} -o build/analyze/anal_v2.o src/analyze/anal_v2.cc
	



