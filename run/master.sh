#! /bin/sh
NRUNS=20
set -e
rm -f -r parameters/run*
rm -f -r parameters/run*
${MADAI_BIN}/parmaker ${NRUNS};
NSTART=0;
qsub -t 1-${NRUNS} ./runpbs.sh ${NSTART}





