Portability
^^^^^^^^^^^

Enabled compiling CUDA device code with clang
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
clang can be used as a device compiler by setting GMX_CLANG_CUDA=ON. A
CUDA toolkit (>=7.0) is also needed. Note that the resulting runtime
performance is usually worse than that of binaries compiled by the
official NVIDIA CUDA compiler (nvcc).

Increased the oldest cmake, compiler and CUDA versions required
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We now require gcc-4.8.1, clang-3.3 and icc-17.0.1, so we can rely on full
C++11 support. We now also require CUDA-6.5 and CMake-3.4.3.

Added check that CUDA available hardware and compiled code are compatible
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added an early check to detect when the :ref:`gmx mdrun` binary does
not embed code compatible with the GPU device it tries to use nor does
it have PTX that could have been just-in-time compiled.

Additionally, if the user manually sets GMX_CUDA_TARGET_COMPUTE=20 and
no later SM or COMPUTE but runs on >2.0 hardware, we'd be executing
just-in-time-compiled Fermi kernels with incorrect host-side code
assumptions (e.g amount of shared memory allocated or texture type).
This change also prevents such cases.

Fixes :issue:`2273`

Disabled ARM Neon native rsqrt iteration used in short-ranged interactions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2261`

Avoided FTZ triggering simd test failures
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
For very small arguments on platforms without FMA support, the Intel
compiler's default usage of flush-to-zero for denormal values can lead
to slight deviations. Since this is a range we really don't care
about, and non-FMA platforms are anyway a thing of the past, just
avoid testing a very small range around that threshold for non-FMA
SIMD platforms.

:issue:`2335`

Fixed OpenCL compiles on Mac OS
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Confirmed to work on Mac OS 10.13.2 running on a Macbook Pro with
Radeon Pro 560.

:issue:`2369`

Tested that nvcc/host compiler combination works
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We now compile a trivial CUDA program during a run of CMake to catch
both unsupported nvcc/host compiler version combinations and other
unknown errors.

:issue:`1616`

Added AVX_512 and KNC symbols to FFTW SIMD test
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Otherwise the CMake code might complain loudly about FFTW not being
accelerated on KNC or KNL hosts.

Implemented changes for CMake policy 0068
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
CMake-3.9 introduced a changed behavior for RPATH vs. install_name
options on OS X. This avoids relying on functionality that will be
removed in future CMake versions.

