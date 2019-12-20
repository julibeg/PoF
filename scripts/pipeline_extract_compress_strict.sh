sbml_fname=${1##*/}; 					#remove dir-part from path
model_name=${sbml_fname%.*}		#remove file extension
BM_rxn=$2

#check if the model has either .xml or .sbml extension
if [ "$extension" != "xml" ] && [ "$extension" != "sbml" ]; then
	echo "ERROR: model file ($1) lacking .sbml or .xml extension"
	exit 1
fi

if [ ! -f $1 ]; then					#check if file exists
	echo "invalid filename: $1"
	exit 2
fi

if [ -d "$model_name" ]; then	#check if directory already exists
	echo "directory with name $model_name already exists"
	exit 3
else
	mkdir $model_name
	mv $1 $model_name
	cd $model_name
fi

printf "$model_name: read model and extract files...\n\n"
extractFromSBML.py $sbml_fname

printf "$model_name: compress model (STDOUT + STDERR in 'compression.log')...\n\n"
compress_network.pl \
	-s ${model_name}.sfile \
	-m ${model_name}.mfile \
	-r ${model_name}.rfile \
	-v ${model_name}.rvfile \
	-p chg_proto.txt \
	-o _comp \
	-l -k \
	> compression.log 2>&1

#get compressed growth reaction and write to tfile_comp
grep -o "[^ ]*$BM_rxn[^ ]*" ${model_name}.rfile_comp \
	> ${model_name}.tfile_comp
printf "$model_name: (compressed) target reaction: "
cat ${model_name}.tfile_comp
echo

echo 	"${model_name}.rfile_comp," \
			"${model_name}.mfile_comp," \
			"${model_name}.sfile_comp," \
			"${model_name}.rvfile_comp," \
			"${model_name}.tfile_comp" > ${model_name}.inp

printf "$model_name: create dual system... (STDOUT + STDERR in 'dual.log')\n\n"
create_ccds_files.pl -c ${model_name}.inp -o ${model_name}_comp \
	> dual.log 2>&1

