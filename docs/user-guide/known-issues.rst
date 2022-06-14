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
