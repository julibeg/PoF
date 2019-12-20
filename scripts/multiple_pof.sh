#!/usr/bin/bash
input_prefix=$1
output_prefix=$2
max_d=$3

for f in ${input_prefix}*; do
	f_num=${f//$input_prefix/}
	#echo -e "${f_num}\t${output_prefix}"
	#PoF.sh $f ${output_prefix}$f_num
	pof=$(PoF.sh $f ${output_prefix}${f_num} $max_d | tail -n3 | head -n1 | awk '{print $3}')
	echo -e "$f_num\t$pof"  
done

