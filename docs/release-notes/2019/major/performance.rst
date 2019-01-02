Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

Implemented update groups
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Domain decomposition can now be based on so-called update groups. These
are groups of atoms with dependencies during the update, which can be
constraints and virtual sites. Update groups can typically be used when
only bonds involving hydrogens are constrained and are enabled
automatically when possible. This improves performance by eliminating
MPI and OpenMP communication for constraints and virtual sites.

PME on GPU when running free energy perturbations not involving charges
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
PME can now be run on a GPU when doing free energy perturbations
that do not involve perturbing charges.

PME long-ranged interaction GPU offload now available with OpenCL
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
On supported devices from all supported vendors (AMD, Intel, NVIDIA),
it is now possible to offload PME tasks to the GPU using OpenCL. This
works in the same way as the former CUDA offload. A single GPU can
now be used to accelerate the computation of the long-ranged PME
interactions. This feature means that only 2-4 CPU cores per
GPU will be about as fast as the 2018 version that needed many more
CPU cores to balance the GPU. Performance on hardware that had good
balance of GPU and CPU also shows minor improvements, and the capacity
for hardware with strong GPUs to run effective simulations is now
greatly improved.

Intel integrated GPUs are now supported for GPU offload with OpenCL
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
On Intel CPUs with integrated GPUs, it is now possible to offload nonbonded tasks
to the GPU the same way as offload is done to other GPU architectures.
This can have performance benefits, in particular on modern desktop and mobile
Intel CPUs this offload can give up to 20% higher simulation performance.

Bonded interactions are now supported for CUDA GPU offload
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Common types of bonded and LJ-14 interactions found can now run on
NVIDIA GPUs with CUDA, with and without domain decomposition.
Interactions with perturbed parameters are not supported.

Added code generation support for NVIDIA Turing GPUs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With CUDA 10.0 NVIDIA Turing GPUs can be directly targeted by the nvcc
compiler. We now generate the appropriate set of flags for the Turing architecture
by default when using CUDA 10 (or later).
