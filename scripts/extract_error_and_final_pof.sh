prefix=$1

echo "d,pof,error,found_sets,all_sets"
for f in ${prefix}*; do
	f_num=${f//$prefix/}
	cutset_line=$(tail -n5 $f | head -n1)
	found_sets=$(echo $cutset_line | awk '{print $5; exit}')
	all_sets=$(echo $cutset_line | awk '{print $6; exit}')
	last_lines=$(tail -n3 $f | head -n4)
	PoF=$(echo $last_lines | awk '{print $3; exit}')
	error=$(echo $last_lines | awk '{print $5; exit}')
	echo "$f_num,$PoF,$error,$found_sets,$all_sets" 
done
