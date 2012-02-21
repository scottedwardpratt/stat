#!/bin/zsh
folders=$(find . -name 'mw*-round' -type d)

for i in $folders; do
		echo "# processing ${i}"
		cd $i
		R --slave < setup-plot-implaus.R
		R --slave < setup-grid-implaus.R
		cd ..
done