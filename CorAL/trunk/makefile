#!/usr/bin/bash
# Could be fixed up in future
export CORAL_HFILES = src/coral/CF2Source/minimization.h\
src/coral/CF2Source/minuitstuff/sfit_minuit.h\
src/coral/CF2Source/sfit.h\
src/coral/CFCalc/cfcalc.h\
src/coral/coral.h\
src/coral/Kernel/kernel.h\
src/coral/Source2CF/source2cf.h\
src/coral/SourceCalc/sourcecalc.h\
src/coral/WaveFunctions/wavefunction.h\
src/coralutils/Arrays/arrays.h\
src/coralutils/GSLMatrix/gslmatrix.h\
src/coralutils/Misc/misc.h\
src/coralutils/NumMath/constants.h\
src/coralutils/NumMath/random.h\
src/coralutils/ParameterMap/parametermap.h\
src/coralutils/SpecialFunctions/Chebyshev/chebyshev.h\
src/coralutils/SpecialFunctions/Chebyshev/chebyshev_fit.h\
src/coralutils/SpecialFunctions/sf.h\
src/coralutils/utilities.h\
src/xgraph/xgraph.h

export CORAL_CCFILES = src/coral/CF2Source/minimization.cc\
src/coral/CF2Source/sfit.cc\
src/coral/CF2Source/sfit_3dgaussian.cc\
src/coral/CF2Source/sfit_blast.cc\
src/coral/CF2Source/sfit_GX1d.cc\
src/coral/CFCalc/cfcalc.cc\
src/coral/Kernel/kernel.cc\
src/coral/Source2CF/source2cf.cc\
src/coral/SourceCalc/sourcecalc.cc\
src/coral/SourceCalc/sourcecalc_blast.cc\
src/coral/SourceCalc/sourcecalc_ellipticblast.cc\
src/coral/SourceCalc/sourcecalc_gauss.cc\
src/coral/SourceCalc/sourcecalc_GX1d.cc\
src/coral/SourceCalc/sourcecalc_OSCAR.cc\
src/coral/WaveFunctions/wf/reid93.cc\
src/coral/WaveFunctions/wf/wf_generic.cc\
src/coral/WaveFunctions/wf/wf_kpluspiminus_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_kpluspiminus_sqwell.cc\
src/coral/WaveFunctions/wf/wf_kpluspiplus_sqwell.cc\
src/coral/WaveFunctions/wf/wf_lambdalambda_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_lambdalambdaantiparspin_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_lambdalambdaparspin_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_nn_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_pipluspiminus_nostrong.cc\
src/coral/WaveFunctions/wf/wf_pipluspiminus_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_pipluspiminus_sqwell.cc\
src/coral/WaveFunctions/wf/wf_pipluspiplus_nostrong.cc\
src/coral/WaveFunctions/wf/wf_pipluspiplus_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_pipluspiplus_sqwell.cc\
src/coral/WaveFunctions/wf/wf_pkplus_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_pkplus_sqwell.cc\
src/coral/WaveFunctions/wf/wf_plambda_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_pn_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_pp_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_pp_schrod.cc\
src/coral/WaveFunctions/wf/wf_ppiminus_sqwell.cc\
src/coral/WaveFunctions/wf/wf_ppiplus_phaseshift.cc\
src/coral/WaveFunctions/wf/wf_ppiplus_sqwell.cc\
src/coral/WaveFunctions/wf/wf_Xipi_phaseshift.cc\
src/coral/WaveFunctions/wfcommon/kpi_phaseshifts.cc\
src/coral/WaveFunctions/wfcommon/partwave.cc\
src/coral/WaveFunctions/wfcommon/pipi_phaseshifts.cc\
src/coral/WaveFunctions/wfcommon/planewave.cc\
src/coral/WaveFunctions/wfcommon/reid93.cc\
src/coral/WaveFunctions/wfcommon/wfcommon.cc\
src/coral/WaveFunctions/wfcommon/wfsquarewell.cc

export CORAL_UTILS_CCFILES = src/coralutils/Arrays/arraycalc.cc\
src/coralutils/Arrays/charray.cc\
src/coralutils/Arrays/mclist.cc\
src/coralutils/Arrays/threed.cc\
src/coralutils/Arrays/ylmarray.cc\
src/coralutils/GSLMatrix/gslmatrix.cc\
src/coralutils/Misc/misc.cc\
src/coralutils/NumMath/random.cc\
src/coralutils/ParameterMap/parametermap.cc\
src/coralutils/SpecialFunctions/Bessel/bess.cc\
src/coralutils/SpecialFunctions/CartesianHarmonics/chcalc.cc\
src/coralutils/SpecialFunctions/CoulWave/coulwave.cc\
src/coralutils/SpecialFunctions/Ylm/ylm.cc

CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}

coral: headers ${MADAI_LIB}/libcoralutils.a ${MADAI_LIB}/libcoral.a ${MADAI_LIB}/libxgraph.a
	for filename in ${CORAL_HFILES};\
	do cp $${filename} ${MADAI_INCLUDE}/ ; done;
	
clean:
	rm ${MADAI_LIB}/libxgraph.a ${MADAI_LIB}/libcoral.a ${MADAI_LIB}/libcoralutils.a

headers: ${CORAL_HFILES}
	for filename in ${CORAL_HFILES};\
	do cp $${filename} include/ ; done;

${MADAI_LIB}/libcoralutils.a : ${CORAL_HFILES} ${CORAL_UTILS_CCFILES}
	mkdir -p build/coralutils;\
	mkdir -p lib;\
	rm *.o build/coralutils/*.o;\
	${CPP} -c -I ${MADAI_GSLPATH}/include -I include -O ${OPT} ${CORAL_UTILS_CCFILES};\
	mv *.o build/coralutils/;\
	ar -ru lib/libcoralutils.a build/coralutils/*.o;\
	mv -f lib/libcoralutils.a ${MADAI_LIB}/

${MADAI_LIB}/libcoral.a : ${CORAL_HFILES} ${CORAL_CCFILES}
	mkdir -p build/coral;\
	mkdir -p lib;\
	rm *.o build/coral/*.o;\
	${CPP} -c -I ${MADAI_GSLPATH}/include -I include -O ${OPT}  ${CORAL_CCFILES};\
	mv *.o build/coral/;\
	ar -ru lib/libcoral.a build/coral/*.o;\
	mv -f lib/libcoral.a ${MADAI_LIB}/

${MADAI_LIB}/libxgraph.a : src/xgraph/xgraph.cc
	mkdir -p build/xgraph;\
	mkdir -p lib;
	rm *.o build/xgraph/*.o;\
	${CPP} -c -I ${MADAI_GSLPATH}/include -I include -I ${MADAI_X11PATH}/include -O ${OPT} src/xgraph/xgraph.cc;\
	mv *.o build/xgraph/;\
	ar -ru lib/libxgraph.a build/xgraph/*.o;\
	mv -f lib/libxgraph.a ${MADAI_LIB}/
