#!/usr/bin/bash
export CORAL_HFILES = minimization.h \
sfit.h \
cfcalc.h \
coral.h \
kernel.h \
source2cf.h \
sourcecalc.h \
wavefunction.h \
arrays.h \
gslmatrix.h \
misc.h \
constants.h \
random.h \
parametermap.h \
chebyshev.h \
chebyshev_fit.h \
sf.h \
utilities.h \
xgraph.h

export CORAL_CCFILES = minimization.cc\
sfit.cc\
sfit_3dgaussian.cc\
sfit_blast.cc\
sfit_GX1d.cc\
cfcalc.cc\
kernel.cc\
kernel_chooser.cc\
source2cf.cc\
sourcecalc.cc\
sourcecalc_blast.cc\
sourcecalc_ellipticblast.cc\
sourcecalc_gauss.cc\
sourcecalc_GX1d.cc\
sourcecalc_OSCAR.cc\
reid93.cc\
wf_generic.cc\
wf_kpluspiminus_phaseshift.cc\
wf_kpluspiminus_sqwell.cc\
wf_kpluspiplus_sqwell.cc\
wf_lambdalambda_phaseshift.cc\
wf_lambdalambdaantiparspin_phaseshift.cc\
wf_lambdalambdaparspin_phaseshift.cc\
wf_nn_phaseshift.cc\
wf_pipluspiminus_nostrong.cc\
wf_pipluspiminus_phaseshift.cc\
wf_pipluspiminus_sqwell.cc\
wf_pipluspiplus_nostrong.cc\
wf_pipluspiplus_phaseshift.cc\
wf_pipluspiplus_sqwell.cc\
wf_pkplus_phaseshift.cc\
wf_pkplus_sqwell.cc\
wf_plambda_phaseshift.cc\
wf_pn_phaseshift.cc\
wf_pp_phaseshift.cc\
wf_pp_schrod.cc\
wf_ppiminus_sqwell.cc\
wf_ppiplus_phaseshift.cc\
wf_ppiplus_sqwell.cc\
wf_Xipi_phaseshift.cc\
kpi_phaseshifts.cc\
partwave.cc\
pipi_phaseshifts.cc\
planewave.cc\
reid93.cc\
wfcommon.cc\
wfsquarewell.cc

export CORAL_UTILS_CC = arraycalc.cc\
charray.cc\
mclist.cc\
threed.cc\
ylmarray.cc\
gslmatrix.cc\
misc.cc\
random.cc\
parametermap.cc\
bess.cc\
chcalc.cc\
chebyshev.cc\
chebyshev_fit.cc\
coulwave.cc\
ylm.cc

CC=/usr/bin/g++
OPT="-fast"

all:
	for filename in ${CORAL_HFILES};\
	do if [ -e include/$${filename} ];\
	then find src -newer include/$${filename} -name $${filename} -exec cp {} include/  \; ;\
	else find src -name $${filename} -exec cp {} include/ \; ;\
	fi;\
	done; \
	for filename in ${CORAL_CCFILES};\
	do if [ -e include/$${filename} ];\
	then find src -newer include/$${filename} -name $${filename} -exec cp {} include/  \; ;\
	else find src -name $${filename} -exec cp {} include/ \; ;\
	fi;\
	done; \
	for filename in ${CORAL_UTILS_CCFILES};\
	do if [ -e include/$${filename} ];\
	then find src -newer include/$${filename} -name $${filename} -exec cp {} include/  \; ;\
	else find src -name $${filename} -exec cp {} include/ \; ;\
	fi;\
	done; \
	
	

	



