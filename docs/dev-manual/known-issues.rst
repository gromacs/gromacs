.. _gmx-dev-known-issues:

Known issues relevant for developers
====================================

This is a non-exhaustive list of known issues that have been observed
and can be of interest for developers. These have not been solved
because they are either outside the scope of the |Gromacs| project
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

PME decomposition automated task assignment broken
--------------------------------------------------

When there are two or more ranks on a node doing combined PP and PME
work (i.e no separate PME ranks) and more GPUs are detected than
ranks, the automated task assignment fails and |Gromacs| aborts with
"Error in user input" message. You can work around
this by using ``-gpu_id`` or ``GMX_GPU_ID`` or limiting the number of
visible GPUs.

:issue:`4684`
