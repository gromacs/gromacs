#!/bin/bash

exeDir="${1:-$(pwd)/../../example}"
export PATH=$exeDir:$PATH

dataDir=${2:-SZ-travis-testdata-master/travis-testdata}
#cp "$exeDir/sz.config" $dataDir
echo "pwd: " $(pwd)
echo "exeDir" ${exeDir}

echo ------------------------------- CESM-ATM-Tylor --------------------------------
app=CESM-ATM-Tylor
errs="1E-1 1E-3 1E-5"
for file in `cd $dataDir/$app;ls *.f32`
do
	filepath=$dataDir/$app/$file
	for err in $errs
	do
		echo sz -z -f -i $filepath -M REL -R $err -2 3600 1800
		sz -z -f -i $filepath -M REL -R $err -2 3600 1800
		echo sz -x -f -i $filepath -s ${filepath}.sz -2 3600 1800 -a
		sz -x -f -i $filepath -s ${filepath}.sz -2 3600 1800 -a
	done
done

echo ------------------------------- EXAFEL --------------------------------
app=EXAFEL
errs="1E-1 1E-3 1E-5"
filepath=$dataDir/$app/smd-cxif5315-r169-calib-fde.f32
for err in $errs
do
	echo sz -z -f -i $filepath -M REL -R $err -3 388 185 320
	sz -z -f -i $filepath -M REL -R $err -3 388 185 320
	echo sz -x -f -i $filepath -s ${filepath}.sz -3 388 185 320 -a
	sz -x -f -i $filepath -s ${filepath}.sz -3 388 185 320 -a
done

cp ../../../example/sz.config .
app=EXAFEL
errs="1E-1 1E-3 1E-5"
filepath=$dataDir/$app/smd-cxif5315-r169-calib-fde.i16
for err in $errs
do
	echo testint_compress -ui16 sz.config $filepath 388 185 320
	testint_compress -ui16 sz.config $filepath 388 185 320
	echo testint_decompress -ui16 ${filepath}.sz 388 185 320
	testint_decompress -ui16 ${filepath}.sz 388 185 320
done

echo ------------------------------- HACC --------------------------------
app=HACC
files="x-131072.f32 y-131072.f32 z-131072.f32"
files2="vx-131072.f32 vy-131072.f32 vz-131072.f32"
errs="1E-1 1E-3 1E-5"

for err in $errs
do
	for file in $files
	do
		filepath=$dataDir/$app/$file
		echo sz -z -f -i $filepath -M REL -R $err -1 131072
		sz -z -f -i $filepath -M REL -R $err -1 131072
		echo sz -x -f -i $filepath -s ${filepath}.sz -1 131072 -a
		sz -x -f -i $filepath -s ${filepath}.sz -1 131072 -a
	done
done

for err in $errs
do
	for file in $files2
	do
		echo sz -z -f -i $filepath -M PW_REL -P $err -1 131072
		sz -z -f -i $filepath -M PW_REL -P $err -1 131072
		echo sz -x -f -i $filepath -s ${filepath}.sz -1 131072 -a
		sz -x -f -i $filepath -s ${filepath}.sz -1 131072 -a
	done
done

echo ------------------------------- Hurricane --------------------------------
app=Hurricane

errs="1E-1 1E-3 1E-5"
for file in `cd $dataDir/$app;ls *.f32`
do
	filepath=$dataDir/$app/$file
	for err in $errs
	do
		echo sz -z -f -i $filepath -M REL -R $err -3 500 500 100
		sz -z -f -i $filepath -M REL -R $err -3 500 500 100
		echo sz -x -f -i $filepath -s ${filepath}.sz -3 500 500 100 -a
		sz -x -f -i $filepath -s ${filepath}.sz -3 500 500 100 -a
	done
done

echo ------------------------------- QMCPACK --------------------------------
app=QMCPack

errs="1E-1 1E-3 1E-5"
for file in `cd $dataDir/$app;ls *.f32`
do
	filepath=$dataDir/$app/$file
	for err in $errs
	do
		echo sz -z -f -i $filepath -M REL -R $err -3 69 69 115 
		sz -z -f -i $filepath -M REL -R $err -3 69 69 115
		echo sz -x -f -i $filepath -s ${filepath}.sz -3 69 69 115 -a
		sz -x -f -i $filepath -s ${filepath}.sz -3 69 69 115 -a
	done
done
