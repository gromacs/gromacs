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
We now require gcc-4.8.1, clang-3.3 and icc-15, so we can rely on full
C++11 support. We now also require CUDA-6.5 and CMake-3.4.3.

Check CUDA available/compiled code compatibility
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

Fix build with cmake 3.10 on gentoo - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

