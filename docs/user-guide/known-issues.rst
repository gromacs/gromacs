Known issues affecting users of |Gromacs|
=========================================

.. _gmx-users-known-issues:

Here is a non-exhaustive list of issues that are we are aware of that are
affecting regular users of |Gromacs|.

Unable to compile with CUDA 11.3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Due to a bug in the nvcc compiler, it is currently not possible
to compile NVIDIA GPU-enabled |Gromacs| with version 11.3 of the CUDA compiler.
We recommend using CUDA 11.4 or newer.

:issue:`4037`

Build is fragile with gcc 7 and CUDA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Different forms of gcc 7 have different behaviour when compiling test
programs with nvcc. This prevents |Gromacs| from reliably testing compilation
flags for use with nvcc. So in this case we use flags unilaterally and this
could lead to compilation errors. The best way to avoid these potential problems
is to use a more recent version of gcc.

:issue:`4478`

SYCL build unstable when using oneAPI with LevelZero backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are multiple issues with different versions of Intel oneAPI when
using the LevelZero backend. 

In many cases, it works fine, and if it fails, it does so explicitly
(either crash or hang), so it should be fine to experiment with.

For most cases, we recommend using OpenCL backend (the default) when
running SYCL build of |Gromacs| on Intel GPUs.

:issue:`4219`
:issue:`4354`

