#!/bin/bash
echo compression

testint_compress -i8 sz_int.config testdata/x86/testint8_8x8x8.dat 8 8 8
testint_decompress -i8 sz_int.config testdata/x86/testint8_8x8x8.dat.sz 8 8 8

testint_compress -i16 sz_int.config testdata/x86/testint16_8x8x8.dat 8 8 8
testint_decompress -i16 sz_int.config testdata/x86/testint16_8x8x8.dat.sz 8 8 8

testint_compress -i32 sz_int.config testdata/x86/testint32_8x8x8.dat 8 8 8
testint_decompress -i32 sz_int.config testdata/x86/testint32_8x8x8.dat.sz 8 8 8

testint_compress -i64 sz_int.config testdata/x86/testint64_8x8x8.dat 8 8 8
testint_decompress -i64 sz_int.config testdata/x86/testint64_8x8x8.dat.sz 8 8 8

testint_compress -ui8 sz_int.config testdata/x86/testint8_8x8x8.dat 8 8 8
testint_decompress -ui8 sz_int.config testdata/x86/testint8_8x8x8.dat.sz 8 8 8

testint_compress -ui16 sz_int.config testdata/x86/testint16_8x8x8.dat 8 8 8
testint_decompress -ui16 sz_int.config testdata/x86/testint16_8x8x8.dat.sz 8 8 8

testint_compress -ui32 sz_int.config testdata/x86/testint32_8x8x8.dat 8 8 8
testint_decompress -ui32 sz_int.config testdata/x86/testint32_8x8x8.dat.sz 8 8 8

