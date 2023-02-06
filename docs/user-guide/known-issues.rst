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

The deform option is not suitable for flow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The deform option currently scales the coordinates, but for flow the deformation
should only be driven by changing periodic vectors. In addition the velocities
of particles need to be corrected when they are displaced by periodic vectors.
Therefore the deform option is currently only suitable for slowly deforming
systems.

:issue:`4607`

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

Unable to build with CUDA 11.5-11.6 and GCC 11 on Ubuntu 22.04
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A bug in the nvcc toolchain, versions 11.5.0-11.6.1, makes it impossible
to build recent |Gromacs| with GCC 11.2 shipped with Ubuntu 22.04. 
We recommend the users to either use an different version of GCC 
(at the time of writing 9.x or 10.x have been reported to work), or manually update the nvcc 
toolchain to version 11.6.2 or newer.

Some non-Ubuntu installations of GCC 11.2 library have been observed to work fine.

When an incompatible combination is used, an error will be raised
from CMake or later during build.

:issue:`4574`

"Cannot find a working standard library" error with ROCm Clang
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some Clang installations don't contain a compatible C++ standard library.
In such cases, you might have to install ``g++`` and help CMake find it
by setting ``-DGMX_GPLUSGPLUS_PATH=/path/to/bin/g++``. 

:issue:`4679`

Ryckaert-Bell dihedral potential calculated imprecisely on Gen9 Intel GPUs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In SYCL/oneAPI build, when bonded forces are offloaded to an Intel Gen9 GPU
(HD Graphics 5xx to 7xx-series; Skylake to Gemini Lake) the Ryckaert-Bell potential
is computed imprecisely. This is unlikely to lead to wrong results, but
we still recommend disabling listed forces offload (``-bonded cpu``) when running
on Gen9 Intel integrated GPUs, especially since offloading is unlikely to offer significant
performance advantage on such devices.

:issue:`4686`

CPU FFT broken with MKL 2023.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When built with MKL 2023.0, CPU FFT fails at plan creation. This is caught by unit tests
(``make check``). We recommend using older version of MKL (2022.3 is known to work)
or :ref:`switching to FFTW3 library <FFT support>`.

:issue:`4691`

Expanded ensemble does not checkpoint correctly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the legacy simulator, because of shortcomings in the
implementation, successful expanded-ensemble MC steps that occured on
checkpoint steps were not recorded in the checkpoint. If that
checkpoint was used for a restart, then it would not necessarily
behave correctly and reproducibly afterwards. So checkpointing of
expanded-ensemble simulations is disabled for the legacy simulator.

Checkpointing of expanded ensemble in the modular simulator works
correctly.

To work around the issue, either avoid ``-update gpu`` (so that it
uses the modular simulator path which does not have
the bug), or use an older version of |Gromacs|
(which does do the buggy checkpointing), or refrain from
restarting from checkpoints in the affected case.

:issue:`4629`
