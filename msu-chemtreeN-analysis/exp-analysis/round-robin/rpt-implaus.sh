#!/bin/zsh
folders=(mw2-round mw3-round mw6-round)

for i in $folders; do
		echo "# processing ${i}"
		cd $i
		R --slave < setup-plot-implaus.R
		R --slave < setup-grid-implaus.R
		cd ..
done