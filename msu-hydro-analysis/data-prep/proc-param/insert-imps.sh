#!/bin/zsh
# insert the impact parameters into the combined params file
AWKFILE=insert-impacts.awk
INFILE=stats.param-combined-s
OUTFILE=params-impact.dat
awk -f ${AWKFILE} < $INFILE > $OUTFILE

