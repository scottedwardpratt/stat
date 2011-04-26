#!/bin/zsh
# ccs
# final step, everything that comes from combine-results is not going to be sorted on ids
# which makes matching against the parameters used to generate the runs a bit more annoying
# so we sort all the data files
echo "sorting  files"
for i in $(find . -name '*combined.dat'); do
		echo $i
		cat $i | sort -n > ${i/.dat}-s.dat
		rm $i
done

for i in $(find . -name '*-errs.dat'); do
		echo $i
		echo ${i/.dat}-s.dat
		cat $i | sort -n > ${i/.dat}-s.dat
		rm $i
done