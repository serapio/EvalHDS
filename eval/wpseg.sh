#! /bin/bash

SEGEVAL=$HOME/Downloads/EvalHDS/texts/segd
testdata=$HOME/Downloads/EvalHDS/wpgold

outHC99dir=$SEGEVAL/HC99_wp
outHCWMdir=$SEGEVAL/HCWM_wp
outHBTdir=$SEGEVAL/HBT_wp
outGBEMdir=$SEGEVAL/GBEM_wp
outBINdir=$SEGEVAL/BIN_wp
outNONEdir=$SEGEVAL/NONE_wp

optdir="wpopts"
mkdir $outHC99dir
mkdir $outHCWMdir
mkdir $outHBTdir $outGBEMdir
mkdir $outBINdir $outNONEdir 
mkdir $optdir

scoresl=wp_walR_eos 
scoresh=wp_wal_eos 
rm $scoresl $scoresh
rm $optdir/*

for file in `ls $testdata` # | head -2`
	do

echo $testdata/$file
echo "filename	$file
gold	$testdata/$file
" >| $optdir/$file

echo "file	$outHC99dir/$file
judge	HC99" >> $optdir/$file

#echo "file	$outHCWMdir/$file
#judge	HCWM" >> $optdir/$file

#echo "file	$outHBTdir/$file
#judge	HBT" >> $optdir/$file

#echo "file	$outGBEMdir/$file
#judge	GBEM" >> $optdir/$file

#echo "file	$outBINdir/$file
#judge	BIN" >> $optdir/$file

echo "file	$outNONEdir/$file
judge	NONE" >> $optdir/$file

	done # end file loop

./Beeferman.py  -wal $optdir/* >> $scoresh
./Beeferman.py  -walR $optdir/* >> $scoresl


#done # end llen loop
