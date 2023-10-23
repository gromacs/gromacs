#!/bin/bash
echo ./szToHDF5 -f sz.config ../../../example/testdata/x86/testfloat_8_8_128.dat 8 8 128
./szToHDF5 -f sz.config ../../../example/testdata/x86/testfloat_8_8_128.dat 8 8 128
echo ./szToHDF5 -d sz.config ../../../example/testdata/x86/testdouble_8_8_128.dat 8 8 128
./szToHDF5 -d sz.config ../../../example/testdata/x86/testdouble_8_8_128.dat 8 8 128
#echo ./szToHDF5 -i8 sz.config ../../../example/testdata/x86/testint8_8x8x8.dat 8 8 8
#./szToHDF5 -i8 sz.config ../../../example/testdata/x86/testint8_8x8x8.dat 8 8 8
#echo ./szToHDF5 -i16 sz.config ../../../example/testdata/x86/testint16_8x8x8.dat 8 8 8
#./szToHDF5 -i16 sz.config ../../../example/testdata/x86/testint16_8x8x8.dat 8 8 8
#echo ./szToHDF5 -i32 sz.config ../../../example/testdata/x86/testint32_8x8x8.dat 8 8 8
#./szToHDF5 -i32 sz.config ../../../example/testdata/x86/testint32_8x8x8.dat 8 8 8
#echo ./szToHDF5 -i64 sz.config ../../../example/testdata/x86/testint64_8x8x8.dat 8 8 8
#./szToHDF5 -i64 sz.config ../../../example/testdata/x86/testint64_8x8x8.dat 8 8 8
echo ./szToHDF5 -u8 sz.config ../../../example/testdata/x86/testint8_8x8x8.dat 8 8 8
./szToHDF5 -u8 sz.config ../../../example/testdata/x86/testint8_8x8x8.dat 8 8 8
echo ./szToHDF5 -u16 sz.config ../../../example/testdata/x86/testint16_8x8x8.dat 8 8 8
./szToHDF5 -u16 sz.config ../../../example/testdata/x86/testint16_8x8x8.dat 8 8 8
echo ./szToHDF5 -u32 sz.config ../../../example/testdata/x86/testint32_8x8x8.dat 8 8 8
./szToHDF5 -u32 sz.config ../../../example/testdata/x86/testint32_8x8x8.dat 8 8 8
echo ./szToHDF5 -u64 sz.config ../../../example/testdata/x86/testint64_8x8x8.dat 8 8 8
./szToHDF5 -u64 sz.config ../../../example/testdata/x86/testint64_8x8x8.dat 8 8 8
