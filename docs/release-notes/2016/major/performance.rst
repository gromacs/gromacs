Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

GPU improvements
^^^^^^^^^^^^^^^^

In addition to those noted below, overall minor improvements contribute
up to 5% increase in CUDA performance, so depending on parameters and compilers
an 5-20% GPU kernel performance increase is expected.
These benefits are seen with CUDA 7.5 (which is now the version we recommend);
certain older versions (e.g. 7.0) see even larger improvements.

Even larger improvements in OpenCL performance on AMD devices are
expected, e.g. can be >50% with RF/plain cut-off and PME with potential shift
with recent AMD OpenCL compilers. 

Note that due to limitations of the NVIDIA OpenCL compiler CUDA is still superior
in performance on NVIDIA GPUs. Hence, it is recommended to use CUDA-based GPU acceleration
on NVIDIA hardware.


Improved support for OpenCL devices
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The OpenCL support is now fully compatible with all intra- and
inter-node parallelization mode, including MPI, thread-MPI, and GPU
sharing by PP ranks. (The previous limitations were caused by bugs in high-level
|Gromacs| code.)

Additionally some prefetching in the short-ranged kernels (similar to
that in the CUDA code) that had been disabled was found to be useful
after all.

Added Lennard-Jones combination-rule kernels for GPUs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Implemented LJ combination-rule parameter lookup in the CUDA and
OpenCL kernels for both geometric and Lorentz-Berthelot combination
rules, and enabled it for plain LJ cut-off. This optimization was
already present in the CPU kernels. This improves performance with
e.g. OPLS, GROMOS and AMBER force fields by about 10-15% (but does not
help with CHARMM force fields because they use force-switched kernels).

Added support for CUDA CC 6.0/6.1
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added build-system and kernel-generator support for the Pascal
architectures announced so far (GP100: 6.0, GP104: 6.1) and supported
by the CUDA 8.0 compiler.

By default we now generate binary as well as PTX code for both sm_60 and
sm_61 and given the considerable differences between the two, we also
generate PTX for both virtual arch. For now we don't add CC 6.2 (GP102)
compilation support as we know nothing about it.

On the kernel-generation side, given the increased register file, for
CC 6.0 the "wider" 128 threads/block kernels are enabled, on 6.1 and
later the 64 threads/block remains.

Improved GPU pair-list splitting to improve performance
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Instead of splitting the GPU lists (to generate more work units) based
on a maximum cut-off, we now generate lists as close to the target
list size as possible. The heuristic estimate for the number of
cluster pairs is now too high by 0-1% instead of 10%. This results in
a few percent fewer pair lists, but still slightly more than
requested.

Improved CUDA GPU memory configuration
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This makes use of the larger amount of L1 cache
available for global load caching on hardware that supports it (K40,
K80, Tegra K1, & CC 5.2) by passing the appropriate command line
option ("-dlcm=ca").

:issue:`1804`

Automatic nstlist changes were tuned for Intel Knight's Landing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

CPU improvements
^^^^^^^^^^^^^^^^

These improvements to individual kernels will provide incremental
improvements to CPU performance for simulations where they are active,
but their value for simulations using GPU offload are much higher,
because via the auto-tuning, they permit all kinds of resource
utilization and throughput to increase.

Optimized the bonded thread force reduction
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The code for multi-threading of bonded interactions has to combine the
forces afterwards. This reduction now uses fixed-size blocks of 32
atoms, and instead of dividing reduction of the whole range of blocks
uniformly over the threads, now only used blocks are divided
(uniformly) over the threads.  This speeds up the reduction by a
factor of the number of threads (!) for typical protein+water systems
when not using domain decomposition. With domain decomposition, the
speed up is up to a factor of 3.

Used SIMD transpose-scatter in bonded force reduction
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The angle and dihedral SIMD functions now use the SIMD transpose
scatter functions for force reduction. This change gives a massive
performance improvement for bondeds, mainly because the dihedral
force update did a lot of vector operations without SIMD that are
now fully replaced by SIMD operations.

Added SIMD implementation of Lennard-Jones 1-4 interactions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The gives a few factors speed improvement. The main improvement comes
from simplified analytical LJ instead of tables; SIMD helps a bit.

Added SIMD implementation of SETTLE
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
On Haswell CPUs, this makes SETTLE a factor 5 faster.

Added SIMD support for routines that do periodic boundary coordinate transformations
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Threading improvements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These improvements enhance the performance of code that runs over
multiple CPU threads.

Improved Verlet-scheme pair-list workload balancing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Implemented near perfect load-balancing for Verlet-scheme CPU
pair-lists. This increases the search cost by 3%, but this is
outweighed by the more balanced non-bonded kernel times, particularly
for small systems.

Improved the threading of virtual-site code
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
On many threads, a significant part of the vsites would end up in
the separate serial task, thereby limiting scaling. Now two weakly
dependent tasks are generated for each thread and one of them uses
a thread-local force buffer, parts of which are reduced by different
threads that are responsible for those parts.

Also the setup now runs multi-threaded.

Add OpenMP support to more loops
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Loops over number of atoms cause significant amount of serial time with
large number of threads, which limits scaling.

Add OpenMP parallelization for the pull code
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The pull code could take up to a third of the compute time for OpenMP
parallel simulation with large pull groups.
Now all pull-code loops over atoms have an OpenMP parallel version.

Other improvements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Multi-simulations are coupled less frequently
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
For example, replica-exchange simulations communicate between simulations
only at exchange attempts. Plain multi-simulations do not communicate
between simulations. Overall performance will tend to improve any time
the progress of one simulation might be faster than others (e.g. it's
at a different pressure, or using a quieter part of the network).
