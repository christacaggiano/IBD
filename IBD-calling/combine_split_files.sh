#!/bin/sh

chr=$1
outlier_file=$2
directory=$3

split -l 1000000 $directory"/chr"$chr".match" $directory"/chr"$chr"_split"

i=0
for f in $directory"/chr"$chr"_split"*;
do
    if [ $i -lt 70 ];
    then

            output_file=$directory"/chr"$chr"_"$i".pkl"

            python combine_pairs_parallel.py $chr $i $f $output_file $outlier_file
    fi
    i=$((i+1))

done
