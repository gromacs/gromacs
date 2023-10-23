#!/bin/bash

if [[ $# < 4 ]]
then
	echo Usage - option 1: $0 [errBoundMode] [error bound] [data directory] [extension] [dimension sizes....]
	echo       - option 2: $0 [errBoundMode] [error bound] [varListFile]
	echo Example: $0 ABS 1E-4 /home/fti/SZ_C_version/CESM-testdata/1800x3600 dat 3600 1800
	exit
fi

errBoundMode=$1
absErrBound=$2

if [ -d $3 ]; then
	option=1
else
	option=0
fi

if [[ $option == 1 ]]; then
	dataDir=$3
	extension=$4
	dim1=$5
	dim2=$6
	dim3=$7
	dim4=$8
else
	varListFile=$3
fi

compressor=sz

#isDimNum is used to indicate the parameter options: either dim1...dim4 are dimensions or dim1 is varList.txt

if [[ $option == 1 ]]; then
	fileList=`cd "$dataDir";ls *.${extension}`
	for file in $fileList
	do
        	echo testfloat_CompDecomp sz.config zc.config "${compressor}($absErrBound)" "$file" $errBoundMode $absErrBound "$dataDir/$file" $dim1 $dim2 $dim3 $dim4
        	./testfloat_CompDecomp sz.config zc.config "${compressor}($absErrBound)" "$file" $errBoundMode $absErrBound "$dataDir/$file" $dim1 $dim2 $dim3 $dim4
	done
else
	nbVars=`./queryVarList -n -i $varListFile`
	for (( i = 0; i < nbVars; i++)); do
		varName=`./queryVarList -m -I $i -i $varListFile`
		file=`./queryVarList -f -I $i -i $varListFile`
		dims=`./queryVarList -d -I $i -i $varListFile`
		echo ./testfloat_CompDecomp sz.config zc.config "${compressor}($absErrBound)" "$varName" $errBoundMode $absErrBound "$file" $dims
		./testfloat_CompDecomp sz.config zc.config "${compressor}($absErrBound)" "$varName" $errBoundMode $absErrBound "$file" $dims
	done
fi

echo "complete"

