#!/bin/bash

mcsfile=$1 	# binary MCS file
outfile=$2
max_d=$3
threads=$4

if [ -z "$max_d" ]
then
	max_d=$(awk -F "" '{ for(i=1; i<=NF;i++) j+=$i; \
		if(j>max) max=j; j=0 } END {print max}' $mcsfile)
fi

if [ -z "$threads" ]
then
	threads=10
fi

failureProbabilityByMcs \
	-i $mcsfile \
	-m $max_d \
	-t $threads \
	-o $outfile
