#!/bin/zsh
folders=(mw2-round mw3-round mw4-round)

for i in $folders; do
		echo "# processing ${i}"
		cd $i
		R --slave < setup.R
		cd ..
done