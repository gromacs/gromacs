Portability
^^^^^^^^^^^

Increased the minimum CUDA version required
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We now require CUDA 7.0, whose features help keep the code more
maintainable.

Increased the minimum MSVC version required
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We now require MSVC 2017, so we can rely on full C++11 support and the
highest quality implementations. On this platform, we now also require
CUDA 9.0.

Updated the OpenCL requirement to version 1.2
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We now require at least OpenCL version 1.2 both for API and kernels. All
currently targeted vendors' libraries do support it, so this is not a
restriction in any way.

Support for ARM Performance Libraries
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The ARM Performance Libraries can now be used for FFT transforms through
the FFTW compatiblity layer, hence GMX_FFT_LIBRARY=fftw3 should be used.
Note that detection is not complete, and therefore
-DFFTWF_LIBRARY=${ARMPL_DIR}/lib/libarmpl_lp64.so
-DFFTWF_INCLUDE_DIR=${ARMPL_DIR}/include
should be manually passed.
