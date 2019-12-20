model_name=$1
dm=$2
d0=$3
t=$4

#set threads for defigueiredo if provided
if [ -z "$t" ]; then
	t=1
fi


printf "$model_name: get MCS... up to d=$dm using $t thread(s) (STDOUT + STDERR in 'defigueiredo.log')\n\n"
defigueiredo \
	-m ${model_name}_comp_dual.mfile \
	-r ${model_name}_comp_dual.rfile \
	-s ${model_name}_comp_dual.sfile \
	-v ${model_name}_comp_dual.vfile \
	-c ${model_name}_comp_dual.cfile \
	-x ${model_name}_comp_dual.xfile \
	-o ${model_name}.mcs.comp \
	-t $t \
	-u $dm \
	-p -i \
	> defigueiredo.log 2>&1

printf "$model_name: transforming compressed MCS to binary representation...\n\n"
mcs2binary.py \
	${model_name}.mcs.comp \
	${model_name}.rfile_comp \
	${model_name}.mcs.comp.binary 

printf "$model_name: calculate PoF up to d=$d0 (STDOUT + STDERR in 'PoF.log')\n\n"
count_compr_rxns.py ${model_name}.rfile_comp ${model_name}.num_comp_rxns
PoFcalc \
	-m ${model_name}.mcs.comp.binary \
	-c ${model_name}.num_comp_rxns \
	-r $(awk '{print NF}' ${model_name}.rfile) \
	-d $d0 \
	-t $t \
	> PoF.log 2>&1

echo "$model_name: done"
