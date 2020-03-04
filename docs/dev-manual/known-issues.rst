.. _gmx-dev-known-issues:

Known issues relevant for developers
====================================

This is a non-exhaustive list of known issues that have been observed
and can be of interest for developers. These have not been solved
because they are either outside the scope of the GROMACS project
or are are simply too difficult or tedious to address ourselves.

Issues with GPU timer with OpenCL
---------------------------------

When building using OpenCL in ``Debug`` mode, it can happen that the GPU timer state gets
corrupted, leading to an assertion failure during the :ref:`mdrun <gmx mdrun>`.
This seems to be related to the load of other, unrelated tasks on the GPU.

GPU emulation does not work
---------------------------

The non-bonded GPU emulation mode does not work, at least for builds
with GPU support; then a GPU setup call is called.
Also dynamic pruning needs to be implemented for GPU emulation.

OpenCL on NVIDIA Volta and later broken
---------------------------------------

The OpenCL code produces incorrect results on Volta and Turing GPU architectures
from NVIDIA (CC 7.0 and 7.5). This is an issue that affects certain flavors of 
the nonboded kernels, most likely a result of miscompilation, and there is no
known workaround.
