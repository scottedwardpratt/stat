#!/bin/zsh

for i in $(find . -name 'mw?-round' -type d); do
		echo "# processing ${i}"
		cd $i
		R --slave < setup.R
		cd ..
done