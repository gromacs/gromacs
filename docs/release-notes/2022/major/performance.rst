Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

GPU direct communication with CUDA-aware MPI
""""""""""""""""""""""""""""""""""""""""""""

Direct GPU communication support has been extended to simulations that use
a CUDA-aware library-MPI when running on NVIDIA GPUs. Detection of CUDA-aware MPI
is performed both at cmake-time and runtime. The feature has been tested
primarily with OpenMPI but any CUDA-aware MPI implementation should be suitable,
and it is also possible to use with the thread-MPI implementation in |Gromacs|.
CUDA-aware MPI support still lacks substantial testing, hence it is included
in the current release as a development feature and should be used with caution.
Hence, even if a suitable MPI is detected, direct communication is not used by
default, but it can be enabled using the GMX_ENABLE_DIRECT_GPU_COMM environment
variable.

:issue:`3960`
:issue:`2915`


Dynamic pairlist generation for energy minimization
"""""""""""""""""""""""""""""""""""""""""""""""""""

With energy minimization, the pairlist, and domain decomposition when running
in parallel, is now performed when at least one atom has moved more than the
half the pairlist buffer size. The pairlist used to be constructed every step.

Nonbonded free-energy kernels use SIMD
""""""""""""""""""""""""""""""""""""""

Free energy calculation performance is improved by making the nonbonded free-energy
kernels SIMD accelerated. On AVX2-256 these kernels are 4 to 8 times as fast.
This should give a noticeable speed-up for most systems, especially if the
perturbed interaction calculations were a bottleneck. This is particularly the
case when using GPUs, where the performance improvement of free-energy runs is
up to a factor of 3.

:issue:`2875`
:issue:`742`

       
PME-PP GPU Direct Communication Pipelining
""""""""""""""""""""""""""""""""""""""""""

For multi-GPU runs with direct PME-PP GPU communication enabled, the
PME rank can now pipeline the coordinate transfers with computation in
the PME Spread and Spline kernel (where the coordinates are
consumed). The data from each transfer is handled separately, allowing
computation and communication to be overlapped. This is expected to
have most benefit on systems where hardware communication interfaces
are shared between multiple GPUs, e.g. PCIe within multi-GPU servers
or Infiniband across multiple nodes.

:issue:`3969`

Domain decomposition with single MPI rank
"""""""""""""""""""""""""""""""""""""""""

When running with a single MPI rank with PME and without GPU, mdrun
will now use the domain decomposition machinery to reorder particles.
This can improve performance, especially for large systems. This
behavior can be controlled with the environment variable
GMX_DD_SINGLE_RANK.

Restricted GPU support with multiple time stepping
""""""""""""""""""""""""""""""""""""""""""""""""""

GPUs can be used in combination with MTS, but for now this is limited
to the setup where only the long-range nonbonded force is applied
in longer timesteps (and computed on the CPU), while all other 
components are are calculated every step (which can be on the GPU).

       
``gmx grompp`` now runs 20-50% faster
"""""""""""""""""""""""""""""""""""""

After a series of improvements, the loops in the parameter- and
atom-lookup code in ``gmx grompp`` have been transformed to
run faster while using simpler, standard code idioms.


PME decomposition support in mixed mode with CUDA and process-MPI
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

PME decomposition is supported now in mixed mode with CUDA backend. 
This is supported only if GROMACS is compiled with external process-MPI 
and underlying MPI implementation is CUDA-aware. This feature lacks substantial testing
and has been disabled by default but can be enabled by setting GMX_GPU_PME_DECOMPOSITION=1 
environment variable.

Performance improvements when running on Ampere-class Nvidia GPUs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Improved performance of the short-ranged non-bonded kernels by up to 12%.

:issue:`3872`
