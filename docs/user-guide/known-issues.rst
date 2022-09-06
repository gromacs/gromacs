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

Verlet buffer underestimated for inhomogeneous systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The current Verlet buffer estimation code assumes that the density
in the system is uniform. This leads to an underestimate of the buffer
for strongly inhomogeneous systems. The temporary solution to this is
to lower the verlet-buffer-tolerance parameter value by the factor between
the uniform density and the local density. In the 2023 release this
correction will be performed automatically.

:issue:`4509`

Verlet buffer underestimated when using only r^-12 potentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When only the repulsive part of the Lennard-Jones potential is used,
as can be the case in coarse-grained systems, the Verlet buffer can be
underestimated due to the extremely non-linear nature of the r^-12 potential.
A temporary solution is to decrease the verlet-buffer-tolerance until you
get a non-zero Verlet buffer. This issue will be fixed in the 2023 release.


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

Unable to build with CUDA 11.6 and gcc-11
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A bug in the nvcc toolchain version 11.6.1 makes it impossible
to build recent |Gromacs| with gcc-11. As these two are the default
versions in Ubuntu 22.04 users are recommended to either install and use
an older version of gcc (version 9.x) has been reported to work, or
manually update the nvcc toolchain to version 11.6.2.

:issue:`4574`

