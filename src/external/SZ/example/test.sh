#If the processor type is in "big endian", please set testdataDir to testdata/ppc instead. See testdata/README for details.
testdataDir=testdata/x86
#make clean
#make

echo ============== testing compression and decompression of 1D array ==============

echo ./testdouble_compress sz.config $testdataDir/testdouble_8_8_128.dat 8192
./testdouble_compress sz.config $testdataDir/testdouble_8_8_128.dat 8192
echo ./testdouble_decompress $testdataDir/testdouble_8_8_128.dat.sz 8192
./testdouble_decompress $testdataDir/testdouble_8_8_128.dat.sz 8192

echo ./testfloat_compress sz.config $testdataDir/testfloat_8_8_128.dat 8192
./testfloat_compress sz.config $testdataDir/testfloat_8_8_128.dat 8192
echo ./testfloat_decompress $testdataDir/testfloat_8_8_128.dat.sz 8192
./testfloat_decompress $testdataDir/testfloat_8_8_128.dat.sz 8192

echo ============== testing compression and decompression of 2D array ==============

echo ./testdouble_compress sz.config $testdataDir/testdouble_8_8_128.dat 64 128
./testdouble_compress sz.config $testdataDir/testdouble_8_8_128.dat 64 128
echo ./testdouble_decompress $testdataDir/testdouble_8_8_128.dat.sz 64 128
./testdouble_decompress $testdataDir/testdouble_8_8_128.dat.sz 64 128

echo ./testfloat_compress sz.config $testdataDir/testfloat_8_8_128.dat 64 128
./testfloat_compress sz.config $testdataDir/testfloat_8_8_128.dat 64 128
echo ./testfloat_decompress $testdataDir/testfloat_8_8_128.dat.sz 64 128
./testfloat_decompress $testdataDir/testfloat_8_8_128.dat.sz 64 128


echo ============== testing compression and decompression of 3D array ==============

echo ./testdouble_compress sz.config $testdataDir/testdouble_8_8_128.dat 8 8 128
./testdouble_compress sz.config $testdataDir/testdouble_8_8_128.dat 8 8 128
echo ./testdouble_decompress $testdataDir/testdouble_8_8_128.dat.sz 8 8 128
./testdouble_decompress $testdataDir/testdouble_8_8_128.dat.sz 8 8 128

echo ./testfloat_compress sz.config $testdataDir/testfloat_8_8_128.dat 8 8 128
./testfloat_compress sz.config $testdataDir/testfloat_8_8_128.dat 8 8 128
echo ./testfloat_decompress $testdataDir/testfloat_8_8_128.dat.sz 8 8 128
./testfloat_decompress $testdataDir/testfloat_8_8_128.dat.sz 8 8 128

echo ./testdouble_compress sz.config $testdataDir/testdouble_8_8_8_128.dat 8 8 8 128
./testdouble_compress sz.config $testdataDir/testdouble_8_8_8_128.dat 8 8 8 128
echo ./testdouble_decompress $testdataDir/testdouble_8_8_8_128.dat.sz 8 8 8 128
./testdouble_decompress $testdataDir/testdouble_8_8_8_128.dat.sz 8 8 8 128

#echo ============== testing batch compression and batch decompression in Fortran =======
#echo ./testdouble_batch_f sz.config
#./testdouble_batch_f sz.config
