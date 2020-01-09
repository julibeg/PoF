#!/bin/bash

mcsfile=$1 	# binary MCS file
outfile=$2
max_d=$3

num_mcs=$(wc -l $mcsfile)
num_rxns=$(awk '{print length; exit}' $mcsfile)

if [ -z "$max_d" ]
then
	max_d=$(awk -F "" '{ for(i=1; i<=NF;i++) j+=$i; \
		if(j>max) max=j; j=0 } END {print max}' $mcsfile)
fi

calcFailureProbability \
	-i $mcsfile \
	-c $num_mcs \
	-r $num_rxns \
	-m $max_d \
	-t 20 \
	-o $outfile
