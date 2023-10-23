SZ2: Error-bounded Lossy Compressor for HPC Data
=====
 (C) 2016-2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
       See COPYRIGHT in top-level directory.

* Major Authors: Sheng Di, Dingwen Tao, Xin Liang, Kai Zhao 
* Supervisor: Franck Cappello 
* Other Contributors: Robert Underwood, Sihuan Li, Ali M. Gok, Cody Rivera, Xiangyu Zou, Wen Xia

## Citations
**Kindly note: This site contains the implementation of SZ2.x. If you mention SZ in your paper, the most appropriate citation is including these three references (***ICDE2021, HPDC2020 and BigData2018***), because they cover the whole design and implementation of the latest version of SZ**.

**Note**: **SZ3** has been released [**here**](https://github.com/szcompressor/SZ3). SZ3 has much higher compression ratios than SZ2 in many cases, with comparable throughput (suffering slightly degraded throughput though). Details can be found in our ICDE21 paper. 

* [SZ3](https://github.com/szcompressor/SZ3): Kai Zhao, Sheng Di, Maxim Dmitriev, Thierry-Laurent D. Tonellot, Zizhong Chen, and Franck Cappello. "[Optimizing Error-Bounded Lossy Compression for Scientiﬁc Data by Dynamic Spline Interpolation](https://ieeexplore.ieee.org/document/9458791)", Proceeding of the 37th IEEE International Conference on Data Engineering (ICDE 21), Chania, Crete, Greece, Apr 19 - 22, 2021.

* SZauto: Kai Zhao, Sheng Di, Xin Liang, Sihuan Li, Dingwen Tao, Zizhong Chen, and Franck Cappello. "[Significantly Improving Lossy Compression for HPC Datasets with Second-Order Prediction and Parameter Optimization](https://dl.acm.org/doi/10.1145/3369583.3392688)", Proceedings of the 29th International Symposium on High-Performance Parallel and Distributed Computing (HPDC 20), Stockholm, Sweden, 2020. (code: https://github.com/szcompressor/SZauto/)

* SZ 2.0+: Xin Liang, Sheng Di, Dingwen Tao, Zizhong Chen, Franck Cappello, "[Error-Controlled Lossy Compression Optimized for High Compression Ratios of Scientific Datasets](https://ieeexplore.ieee.org/document/8622520)", in IEEE International Conference on Big Data (Bigdata 2018), Seattle, WA, USA, 2018.

* SZ 1.4.0-1.4.13: Dingwen Tao, Sheng Di, Franck Cappello, "[Significantly Improving Lossy Compression for Scientific Data Sets Based on Multidimensional Prediction and Error-Controlled Quantization](https://ieeexplore.ieee.org/document/7967203)", in IEEE International Parallel and Distributed Processing Symposium (IPDPS 2017), Orlando, Florida, USA, 2017.

* SZ 0.1-1.0: Sheng Di, Franck Cappello, "[Fast Error-bounded Lossy HPC Data Compression with SZ](https://ieeexplore.ieee.org/document/7516069)", in IEEE International Parallel and Distributed Processing Symposium (IPDPS 2016), Chicago, IL, USA, 2016.

* Point-wise relative error bound mode (i.e., PW_REL): Xin Liang, Sheng Di, Dingwen Tao, Zizhong Chen, Franck Cappello, "[An Efficient Transformation Scheme for Lossy Data Compression with Point-wise Relative Error Bound](https://ieeexplore.ieee.org/document/8514879)", in IEEE International Conference on Clustering Computing (CLUSTER 2018), Belfast, UK, 2018. (Best Paper)

This document simply introduces how to install and use the SZ compressor. More details can be found in doc/userguide.pdf. 

*OpenCL version can be found in the package, while this is a deprecated code for **GPU**. The optimized GPU code in CUDA can be found at [https://github.com/szcompressor/cuSZ](https://github.com/szcompressor/cuSZ).*

## Installation

### Installation way 1:
* ./configure --prefix=[INSTALL_DIR] (Please use --enable-fortran if you need Fortran interface)
* make
* make install

### Installation way 2:
* mkdir build && cd build
* cmake .. -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR]
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and .a and .so libraries in [INSTALL_DIR]/lib

## Testing Examples
--------------------------------------

Examples can be found in the [SZ_PACKAGE]/example

You can use the executable 'sz' command to do the compression/decompression. Please see the user guide or run 'sz --help' for details.

Alternatively, you can also also call our API to do the compression/decompressoin. Here are two examples: testfloat_compress.c and testfloat_decompress.c

## Compression
--------------
* ./test_compress sz.config testdouble_8_8_8_128.dat 8 8 8 128
* ./test_compress sz.config testdouble_8_8_128.dat 8 8 128

`Decription: `

testdouble_8_8_128.dat and testdouble_8_8_8_128.dat are two binary testing files (small-endian format), which contains a 3d array (128X8X8) and a 4d array (128X8X8X8) respectively. Their data values are shown in the two plain text files, testdouble_8_8_128.txt and testdouble_8_8_8_128.txt. These two data files are from FLASH_Blast2 and FLASH_MacLaurin respectively (both are at time step 100). The compressed data files are namely testdouble_8_8_8_128.dat.sz and testdouble_8_8_128.dat.sz respectively.

sz.config is the configuration file. The key settings are errorBoundMode, absErrBound, and relBoundRatio, which are described below.

* absErrBound refers to the absolute error bound, which is to limit the (de)compression errors to be within an absolute error. For example, absErrBound=0.0001 means the decompressed value must be in [V-0.0001,V+0.0001], where V is the original true value.

* relBoundRatio refers to relative bound ratio, which is to limit the (de)compression errors by considering the global data value range size (i.e., taking into account the range size (max_value - min_value)). For example, suppose relBoundRatio is set to 0.01, and the data set is {100,101,102,103,104,...,110}, so the global value range size is 110-100=10, so the error bound will actually be 10*0.01=0.1, from the perspective of "relBoundRatio".

* errorBoundMode is to indicate the error-bounding way in the compression, such as based on absolute error bound, relative error bound, etc. 
The options are shown below.
	* ABS: take only "absolute error bound" into account. That is, relative bound ratio will be ignored.
	* REL: take only "relative bound ratio" into account. That is, absolute error bound will be ignored. 
	* ABS_AND_REL: take both of the two bounds into account. The compression errors will be limited using both absErrBound and relBoundRatio*rangesize. That is, the two bounds must be both met.
	* ABS_OR_REL: take both of the two bounds into account. The compression errors will be limited using either absErrBound or relBoundRatio*rangesize. That is, only one bound is required to be met.
	* PW_REL: take "point-wise relative error bound" in the compression. 

## Decompression

* ./test_decompress testdouble_8_8_8_128.dat.sz 8 8 8 128
* ./test_decompress testdouble_8_8_128.dat.sz 8 8 128

The output files are testdouble_8_8_8_128.dat.sz.out and testdouble_8_8_128.dat.sz.out respectively. You can compare .txt file and .out file for checking the compression errors for each data point. For instance, compare testdouble_8_8_8_128.txt and testdouble_8_8_8_128.dat.sz.out.

## Application Programming Interface (API)

Programming interfaces are procides in two programming languages - C: 

The interfaces are listed below. More details can be found in the user guide. 

### C Interface

* int SZ_Init();

* char *SZ_compress(void *data, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
* char *SZ_compress_args(void *data, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

* double *SZ_decompress(char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
* void SZ_Finalize();

### Fortran Interface

Please see doc/use-guide for details

### Python Interface

**NOTE: THESE BINDINGS ARE DEPRECATED **

The following information is provided for historical purposes only.
Please consider updating to using the Python bindings for SZ provided with [LibPressio](https://github.com/codarcode/libpressio#python) instead which are more efficient and updated with new features in SZ as they are developed.

The python bindings requires some additional dependencies:

- python with development libraries
- numpy with development libraries
- swig-3.0.12 or newer
- cmake-3.14 or newer

To use the python interface, autotools is not supported.  You should compile with CMake instead as follows:

```bash
# if you are working in a python virtual envionment, source it here
# otherwise cmake will detect the wrong numpy version which can cause
# segmentation faults and other bizzare errors
source bin/actiate

mkdir build
cd build
#you can specify other cmake arguments such as CMAKE_INSTALL_PREFIX here
cmake .. -DBUILD_PYTHON_WRAPPER=ON
cmake --build .

#for system wide installation sudo is required
sudo cmake --install .
```

Please note, building static libraries is incompatiable with building the python wrappers, and is building the python wrappers is disabled if `-DBUILD_SHARED_LIBS=OFF`

An example usage file can be found in `example/test.py`

Additional documentation can be found using the `help` function in python.

## Limitation of this version

SZ is not suitable for compressing tiny datasets (such as the size <10KB)

## GPU Version

Please refer to [this repository](https://github.com/szcompressor/cuSZ) for our GPU/CUDA version of SZ (called cuSZ). Please create an issue ticket there if you have any questions or issues regarding cuSZ.

## Performance Portable Version

Please refer to [this repository](https://github.com/szcompressor/kokkosSZ) for our performance portable version of SZ (called kSZ) using [Kokkos](https://github.com/kokkos/kokkos) programming model. Please create an issue ticket there if you have any questions or issues regarding kSZ.

## Version history

version		New features

* SZ 0.2-0.4	Compression ratio is the same as SZ 0.5. The key difference is different implementation ways, such that SZ 0.5 is much faster than SZ 0.2-0.4.
* SZ 0.5.1	Support version checking
* SZ 0.5.2	finer compression granularity for unpredictable data, and also remove redundant Java storage bytes
* SZ 0.5.3 	Integrate with the dynamic segmentation support
* SZ 0.5.4	Gzip_mode: defaut --> fast_mode ; Support reserved value	
* SZ 0.5.5	runtime memory is shrinked (by changing int xxx to byte xxx in the codes). The bug that writing decompressed data may encounter exceptions is fixed. Memory leaking bug for ppc architecture is fixed.
* SZ 0.5.6	improve compression ratio for some cases (when the values in some segementation are always the same, this segment will be merged forward)
* SZ 0.5.7	improve the decompression speed for some cases
* SZ 0.5.8	Refine the leading-zero granularity (change it from byte to bits based on the distribution). For example, in SZ0.5.7, the leading-zero is always in bytes, 0, 1, 2, or 3. In SZ0.5.8, the leading-zero part could be xxxx xxxx xx xx xx xx xxxx xxxx (where each x means a bit in the leading-zero part)
* SZ 0.5.9	optimize the offset by using simple right-shifting method. Experiments show that this cannot improve compression ratio actually, because simple right-shifting actually make each data be multiplied by 2^ {-k}, where k is # right-shifting bits. The pros is to save bits because of more leading-zero bytes, but the cons is much more required bits to save. A good solution is SZ 0.5.10!
* SZ 0.5.10	optimze the offset by using the optimized formula of computing the median_value based on optimized right-shifting method. Anyway, SZ0.5.10 improves compression ratio a lot for hard-to-compress datasets. (Hard-to-compress datasets refer to the cases whose compression ratios are usually very limited)
* SZ 0.5.11	In a very few cases, SZ 0.5.10 cannot guarantee the error-bounds to a certain user-specified level. For example, when absolute error bound = 1E-6, the maximum decompression error may be 0.01(>>1E-6) because of the huge value range even in the optimized segments such that the normalized data cannot reach the required precision even stoaring all of the 64 or 32 mantissa bits. SZ 0.5.11 fixed the problem well, with degraded compression ratio less than 1%.
* SZ 0.5.12 	A parameter setting called "offset" is added to the configuration file sz.config. The value of offset is an integer in [1,7]. Generally, we recommend offset=2 or 3, while we also find that some other settings (such as offset=7) may lead to better compression ratios in some cases. How to automize/optimize the selection of offset value would be the future work. In addition, the compression speed is improved, by replacing java List by array implementation in the code.
* SZ 0.5.13	Compression performance is improved, by replacing some class instances in the source code by primitive data type implementation.
* SZ 0.5.14 	fixed a design bug, which improves the compression ratio further.
* SZ 0.5.15	improved the compression ratio for single-precision data compression, by tuning the offset.

The version 0.x were all coded in Java, and C/Fortran interfaces were provided by using JNI and C/Fortran wrapper. SZ 1.0 is coded in C purely.

* SZ 1.0		Pure C version. In this version, the users don't need to install JDK and make the relative configurations any more. It provides dataEndienType in the sz.config 		file, so it can be used to compress the data file which was generated on different endian-type systems.
* SZ 1.2		Improve compression: using 8 intervals (best selection between linear-curve-fitting and quadratic-curve fitting)
* SZ 1.3		Significantly improved compression: using 255 intervals, with multi-dimensional data prediction (on 1D, 2D, and 3D)
* SZ 1.4		Use 65536 intervals
* SZ 1.4.2	Extending the number of intervals from 255 to 65536, by tailoring/reimplementing the Huffman encoding by ourselves.
* SZ 1.4.3	Add the intervals_count to the configuration file (sz.config), allowing users to control it.
* SZ 1.4.4	Remove segmentation step quantization_intervals
* SZ 1.4.5	Optimize the number of intervals: the # intervals will be automatically optimized before the compression if quantization_itnervals is set to 0.
* SZ 1.4.6-beta   Three compression modes are provided (SZ_BEST_SPEED, SZ_BEST_COMPRESSION, SZ_DEFAULT_COMPRESSION), the maximum # quantization intervals is 65536.
* SZ 1.4.7-beta   Fix some mem leakage bugs. Fix the bugs about memory crash or segmentation faults when the number of data points is pretty large. Fix the sementation fault bug happening when the data size is super small for 1D array. Fix the error bound may not be guaranteed in some cases.
* SZ 1.4.9-beta	Support point-wise relative error bound setting, and optional Fortran compilation.
* SZ 1.4.9.1-beta Fix the bug in the fortran interface about SZ_batch_compression
* SZ 1.4.9.2-beta Update the user guide by describing how to optimize the compression quality on demand.
* SZ 1.4.9.3-beta Fix the sementation fault bug happening when the data size is super small for 2D array and 3D array. (Specifically, when the data size is very small while the error bound is set to very small too, the Huffman tree overhead will be relatively huge such that the compressed size may exceed the original data size, leading to segmentation fault when further compressing it by the last lossless compression step. Solution: In this case, the data will be compressed by Zlib for simplicity, with no compression errors).
* SZ 1.4.10-beta (1) Support direct sub-block data compression; (2) Support compression of large data file directly (i.e., the number of data points could be up to as large as LONG size, unlike the previous version that can only compress 2^{32} data points each time): that is, int nbEle --> size_t nbEle; (3) separate the internel functions from the sz.h; 
* SZ 1.4.11-beta (1) Support HDF5. (2) Support integer data compression. (3) Provide optional wavelet transform as a preprocessing step in SZ and an optional Tucker tensor decomposition.
* SZ 1.4.11: (1) This is a stable version which have went through a long period test. (2) Fix a small bug (the maximum compression error may be slightly greater than error bound in some cases); (3) Support integer compression (for all types of integers); (4) Support HDF5-SZ for all types of integers; (5) Support getting the metadata from a given compressed data file (by using SZ_getMetaData) and printing the metadata (by SZ_printMetaData) with -p option of the executable command ("sz"); (6) Change libsz.a to libSZ.a in case of conflict with szip (note that szip has already used libsz.a).
* SZ 1.4.12: (1) Support thresholding-based strategy for 1D data compression based on point-wise relative error bound. (In order to test it, please select errBoundMode = PW_REL, and set the point-wise relative error bound using the parameter pw_relBoundRatio in the sz.config.) For other dimensions of data, point-wise relative error based compression is using block-based strategy (see our DRBSD-2 paper for details) (2) fix the bug in the callZlib.c (previously, segmentation fault might happen when using best_compression mode). (2) Fix a small bug that happened when the data size is extremely huge (nbEle>4G) and the compression mode is SZ_BEST_COMPRSSION. Specifically, the previous call to zlib functions has one potential bug that may lead to segmentation fault, which has been fixed.
* SZ 1.4.12.1: (1) Fix the bug tha segmentation fault may happen when the error bound is greater than the value range size. 
* SZ 1.4.13.0: Support openMP
* SZ 1.4.13.1: Integrate PaSTRI algorithm into SZ for GAMESS two-electron integral data compression.
* SZ 1.4.13.2: Clean the code by putting all the global variables into two data structures (conf_params and exe_params) and change Huffman tree encoder from global state to local state in order to support thread-safe feature.
* SZ 1.4.13.3: Support the time-based compression for 1D array
* SZ 1.4.13.4: Support the time-based compression for 1D,2D, and 3D array (use --enable-timecmpr to switch it on). SZ 1.4.13.0/1/2/3/ has a bug when doing compression/decompression alternatively in the same progress. SZ 1.4.13.4 fixed it.
* SZ 1.4.13.5: remove the pwrType in SZ_compress_args() interface, because of uselessness of pwrType setting in the new design of point-wise relative error based compression in the later version of SZ.
* SZ 2.0.0.0: Significantly improve the compression quality for the high-compression cases (i.e., improve the PSNR for the cases with high compression ratios)
* SZ 2.0.1.0: Further improve the compression quality for high-compression cases than 2.0.0.0. Moreover, improve point-wise relative error bounded compression.
* SZ 2.0.2.0: Further improve the compression/decompression rate and also compression ratio (by 10-20%), by replacing Zlib by Zstd as default setting.
* SZ 2.1: Significantly improve the compressoin speed for point-wise relative error bound (please see our paper published in MSST19 for details)
* SZ 2.1.1: We accelerated compression for absolute error bound slightly. We fixed the compilation issue for OSX.
* SZ 2.1.2: Fix memory a little memory leak for point-wise relative error bound Remove a useless step (optimization of # quantization bins) for point-wise relative error bound mode, such that the compression rate is further improved (3%).
* SZ 2.1.4.1: Fix a few bugs which were introduced when integrating fast point-wise relative error bounded compression. The bugs could cause unbounded decompression and slow speed, especially for 1D datasets. CMakeLists are also revised to enable SZ to be compiled on Cray system such as Cori.
* SZ 2.1.4.2: Revise a bug in 2.1.4.2: The bug happened only when the data size is too small (smaller than a datablock such as 6x6x6) for 2D and 3D cases.
* SZ 2.1.5: Fix some bugs in temporal compression and support random-access variable selection for compression/decompression. Fix compilation bug for Fortran version. 
* SZ 2.1.6: Fix a bug (error couldn't be bounded) happening when setting point-wise relative error bound for 3D double-precision data.
