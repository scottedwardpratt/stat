#! /bin/sh
NRUNS=5
set -e
rm -f -r parameters/run*
${MADAI_BIN}/parmaker ${NRUNS}
irun=1;
while [ ${irun} -le ${NRUNS} ]
do
	echo irun=${irun};
	./isHydro3 run${irun};
	./b3d run${irun};
	irun=`expr ${irun} + 1`;
done




