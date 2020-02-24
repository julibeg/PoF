#!/bin/bash

sbml_fname=${1##*/}; 						#remove dir-part from path
extension="${sbml_fname##*.}"		#store extension
model_name=${sbml_fname%.*}			#remove file extension
BM_rxn=$2
dm=$3
d0=$4
t=1

#check if the model has either .xml or .sbml extension
if [ "$extension" != "xml" ] && [ "$extension" != "sbml" ] && [ "$extension" != "json" ]; then
	echo "ERROR: model file ($1) lacking .sbml/.xml/.json extension"
	exit 1
fi

#set threads for defigueiredo if provided
if [ ! -z "$5" ]; then
	t=$5
fi

if [ ! -f $1 ]; then							#check if file exists
	echo "ERROR: invalid filename: $1"
	exit 2
fi

if [ -d "$model_name" ]; then			#check if directory already exists
	echo "ERROR: directory with name $model_name already exists"
	exit 3
else
	mkdir $model_name
	mv $1 $model_name
	cd $model_name
fi

printf "$model_name: read model and extract files...\n\n"
extractFromSBML.py $sbml_fname

#biomass rxn won't be compressed
printf "$model_name: biomass reaction is $BM_rxn\n\n"
echo "$BM_rxn" > ${model_name}.nfile

printf "$model_name: compress model (STDOUT + STDERR in 'compression.log')...\n\n"
(time compress_network.pl \
	-s ${model_name}.sfile \
	-m ${model_name}.mfile \
	-r ${model_name}.rfile \
	-v ${model_name}.rvfile \
	-i ${model_name}.nfile \
	-p chg_proto.txt \
	-o _comp \
	-l -k) \
	> compression.log 2>&1

#write growth rxn to tfile_comp
cp ${model_name}.nfile ${model_name}.tfile_comp


echo 	"${model_name}.rfile_comp," \
			"${model_name}.mfile_comp," \
			"${model_name}.sfile_comp," \
			"${model_name}.rvfile_comp," \
			"${model_name}.tfile_comp" > ${model_name}.inp

printf "$model_name: create dual system... (STDOUT + STDERR in 'dual.log')\n\n"
(time create_ccds_files.pl -c ${model_name}.inp -o ${model_name}_comp) \
	> dual.log 2>&1

printf "$model_name: get MCS... up to d=$dm using $t thread(s) (STDOUT + STDERR in 'defigueiredo.log')\n\n"
(time defigueiredo \
	-m ${model_name}_comp_dual.mfile \
	-r ${model_name}_comp_dual.rfile \
	-s ${model_name}_comp_dual.sfile \
	-v ${model_name}_comp_dual.vfile \
	-c ${model_name}_comp_dual.cfile \
	-x ${model_name}_comp_dual.xfile \
	-o ${model_name}.mcs.comp \
	-t $t \
	-u $dm \
	-p -i) \
	> defigueiredo.log 2>&1

printf "$model_name: transforming compressed MCS to binary representation...\n\n"
mcs2binary.py \
	${model_name}.mcs.comp \
	${model_name}.rfile_comp \
	${model_name}.mcs.comp.binary 

printf "$model_name: calculate PoF up to d=$d0 (STDOUT + STDERR in 'PoF.log')\n\n"
count_compr_rxns.py ${model_name}.rfile_comp ${model_name}.num_comp_rxns
(time PoFcalc \
	-m ${model_name}.mcs.comp.binary \
	-c ${model_name}.num_comp_rxns \
	-r $(awk '{print NF}' ${model_name}.rfile) \
	-d $d0 \
	-t $t) \
	> PoF.log 2>&1

echo "$model_name: done"
