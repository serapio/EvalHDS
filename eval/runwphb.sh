#! /bin/bash

# This script demonstrates the procedure for segmenting the WP corpus with HierBayes

for file in `ls ../wpeos_hb/*.txt.eos.tok`
do 
	fname=`basename $file`
	echo $fname
	./hsegment 3 hconfig/bayes-lsvb-top.config < ../wpeos_hb/$fname > ../wpeos_hb-top/$fname
	#./hsegment hconfig/greedy-bayes.config < ../wpeos_hb/$fname > ../wpeos_gb/$fname
done

python ./hb-out-convert.py ../wpeos_hb-top ../wpeos_noseg ../wpeos_hb-out
#python ./hb-out-convert.py ../wpeos_gb ../wpeos_noseg ../wpeos_gb-out

