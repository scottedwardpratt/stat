#! /bin/sh
NRUNS=10
set -e
rm -f -r parameters/run*
rm -f -r output/run*
${MADAI_BIN}/parmaker ${NRUNS};
qsub -t 1-${NRUNS} runpbs.sh





