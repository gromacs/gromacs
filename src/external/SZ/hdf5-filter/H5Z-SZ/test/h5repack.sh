#!/bin/bash

if [[ $# != 2 ]]
then
	echo "Usage: $0 [inputf_hdf5_ile] [output_hdf5_file]"
	echo "Example: $0 testfloat_8_8_128.h5 testfloat_8_8_128_sz.h5"
	exit
fi

inputFile=$1
outputFile=$2
h5repack -f UD=32017,0 -i $inputFile -o $outputFile
