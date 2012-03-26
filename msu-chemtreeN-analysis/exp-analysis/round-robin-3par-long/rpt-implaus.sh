#!/bin/zsh

for i in $(find . -name 'mw*-round' -type d); do
		echo "# processing ${i}"
		cd $i
		R --slave < setup-plot-implaus.R
		R --slave < setup-grid-implaus.R
		cd ..
done