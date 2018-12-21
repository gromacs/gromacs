.. _gmx-dev-known-issues:

Known issues relevant for developers
====================================

This is a non-exhaustive list of known issues that have been observed
and can be of interest for developers. These have not been solved
because they are either outside the scope of the GROMACS project
or are are simply too difficult or tedious to address ourselves.

FP exceptions with CUDA 7.0
---------------------------

When using CUDA 7.0 in a ``Debug`` build, if the PME FFT task is offloaded
to a GPU, a floating point exception will abort the :ref:`mdrun <gmx mdrun>` execution.
The exception originates from the CUDA FFT (cuFFT) library.
To avoid this issue, we advise using a later CUDA version.

Issues with GPU timer with OpenCL
---------------------------------

When building using OpenCL in ``Debug`` mode, it can happen that the GPU timer state gets
corrupted, leading to an assertion failure during the :ref:`mdrun <gmx mdrun>`.
This seems to be related to the load of other, unrelated tasks on the GPU.

