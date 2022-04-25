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

