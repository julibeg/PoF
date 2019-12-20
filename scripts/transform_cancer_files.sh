#!/bin/bash

model_name=$1

cat ${model_name}.Reactions | tr '\n' ' ' > ${model_name}.rfile
echo >> ${model_name}.rfile

cat ${model_name}.Metabolites | tr '\n' ' ' > ${model_name}.mfile
echo >> ${model_name}.mfile

cat ${model_name}.rev | tr '\n' ' ' > ${model_name}.rvfile
echo >> ${model_name}.rvfile

cp ${model_name}.S ${model_name}.sfile

grep -i "biomass" ${model_name}.Reactions | tr '\n' ' ' > ${model_name}.nfile
echo >> ${model_name}.nfile
