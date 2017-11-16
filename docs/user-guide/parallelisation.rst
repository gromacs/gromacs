.. _gmx-parallel:

Acceleration and parallelization
================================

Here we give an overview on the parallelization and acceleration schemes employed by |Gromacs|.
The aim is to provide an understanding of the underlying mechanism that make |Gromacs| one of the
fastest molecular dynamics packages. On the other hand, the information presented
should help choosing the appropriate parallelization options, run configuration,
as well as acceleration options in order to achieve optimal simulation performance.

Acceleration
------------

CPU acceleration: SSE, AVX, etc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In |Gromacs|, three levels of non-bonded kernels are used: reference or generic
kernels, optimized plain-C kernels and SIMD intrinsics accelerated kernels. Other compute
intensive parts of the code, mainly PME, bonded force calculation, and neighbour searching
also employ SIMD intrinsic acceleration.

The SIMD intrinsic code is compiled by the compiler.
Technically, it is possible to compile different levels of acceleration into one binary,
but this is difficult to manage with acceleration in many parts of the code. Thus, you need
to configure and compile |Gromacs| with a single target hardware acceleration which corresponds
to a SIMD instruction set. By default, the build system will detect the highest supported
acceleration of the host where the compilation is carried out. For cross-compiling for 
a machine with a different highest SIMD instructions set, in order to set the target acceleration,
the ``-DGMX_SIMD`` CMake option can be used. For best performance always pick the highest
(latest) SIMD instruction set supported by the target architecture (and |Gromacs|). To use a single
installation on multiple different machines, it is convenient to compile the analysis tools with
the lowest common SIMD instruction set (as these rely little on SIMD acceleration), but for best
performance :ref:`mdrun <gmx mdrun>` should be compiled separately for each machine.

Additional information can be found in the installation guide under :ref:`SIMD support <gmx-simd-support>`.

GPU acceleration
^^^^^^^^^^^^^^^^

The most compute-intensive parts of simulations, non-bonded force calculation, as well
as possibly the PME and bonded force calculation can be
offloaded a GPU and carried out simultaneously with CPU calculations for constraints and
position updates. Native GPU acceleration is supported with the
:doc:`Verlet cut-off scheme <cutoff-schemes>`
(not with the group scheme) with PME, reaction-field, and plain cut-off electrostatics.    

For more information about the GPU kernels, please see the :ref:`Installation guide <gmx-gpu-support>`.

The native GPU acceleration can be turned on or off, either at run-time using the
:ref:`mdrun <gmx mdrun>` ``-nb`` option, or at configuration time using the ``GMX_GPU`` CMake variable.
 
Parallelization schemes
-----------------------

|Gromacs|, being performance-oriented, has a strong focus on efficient parallelization.
There are multiple parallelization schemes available, therefore a simulation can be run on a
given hardware with different choices of run configuration.

MPI
^^^

Parallelization based on MPI has been part of |Gromacs| from the early versions hence is compatible
with the majority of MD algorithms. At the heart of the MPI parallelization is the neutral-territory
:ref:`domain decomposition <gmx-domain-decomp>` which supports fully automatic dynamic load balancing.

To parallelize simulations across multiple machines (e.g. nodes of a cluster) 
:ref:`mdrun <gmx mdrun>` needs to be compiled with MPI which can be enabled using the ``GMX_MPI`` CMake variable.
 
Multithreading with thread-MPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The thread-MPI library provides an implementation of a subset of the MPI 1.1 specification,
based on the system threading support. Both POSIX pthreads and Windows threads are supported,
thus the implementation provides great portability to most UNIX/Linux and Windows operating systems.
Acting as a drop-in replacement for MPI, thread-MPI enables compiling and running :ref:`mdrun <gmx mdrun>`
on a single machine (i.e. not across a network) without MPI. Additionally, it not only provides a
convenient way to use computers with multicore CPU(s), but thread-MPI does in some
cases make :ref:`mdrun <gmx mdrun>` run slightly faster than with MPI.
 
Thread-MPI is included in the |Gromacs| source and it is the default parallelization since
version 4.5, practically rendering the serial :ref:`mdrun <gmx mdrun>` deprecated.
Compilation with thread-MPI is controlled by the ``GMX_THREAD_MPI`` CMake variable.
 
Thread-MPI is compatible with most :ref:`mdrun <gmx mdrun>` features and parallelization schemes, 
including OpenMP, GPUs; it is not compatible with MPI and multi-simulation runs.
 
By default, the thread-MPI mdrun will use all available cores in the machine by starting
as many ranks as the number of cores. The number of ranks can be controlled using the
``-nt`` and ``-ntmpi`` options. ``-nt`` represents the total number of threads
to be used (which can be a mix of thread-MPI and OpenMP threads with the
:doc:`Verlet scheme <cutoff-schemes>`).
 
Multi-level parallelization: MPI and OpenMP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The multi-core trend in CPU development substantiates the need for multi-level parallelization.
Current multiprocessor machines can have 2-4 CPUs with a core count as high as 64. As the memory
and cache subsystem is lagging more and more behind the multicore evolution, this emphasizes
non-uniform memory access (NUMA) effects, which can become a performance bottleneck. At the same
time, all cores share a network interface. In a purely MPI-parallel scheme, all MPI processes
use the same network interface, and although MPI intra-node communication is generally efficient,
communication between nodes can become a limiting factor to parallelization. This is especially
pronounced in the case of highly parallel simulations with PME (which is very communication
intensive) and with ``fat`` nodes connected by a slow network. Multi-level parallelism aims
to address the NUMA and communication related issues by employing efficient
intra-node parallelism, typically multithreading.
 
OpenMP multithreading is supported in :ref:`mdrun <gmx mdrun>`
and combined with MPI (or thread-MPI) it enables multi-level and heterogeneous parallelization.
With the :doc:`Verlet cut-off scheme <cutoff-schemes>` full OpenMP multithreading support is implemented,
but the group scheme currently only supports OpenMP threading for PME. 
 
OpenMP is enabled by default in |Gromacs| and can be turned on/off at configure time with
the ``GMX_OPENMP`` CMake variable and at run-time with the ``-ntomp`` option (or
``OMP_NUM_THREADS`` enviroenment variable).
 
While the OpenMP implementation itself is quite efficient and scales well (up to 12-24 threads
on Intel and 6-8 threads on AMD CPUs), when combining with MPI it has an additional overhead
especially when running separate multi-threaded PME nodes. Depending on the architecture,
input system size, as well as other factors, MPI+OpenMP runs can be as fast and faster
already at small number of processes (e.g. multi-processor Intel Westmere or Sandy Bridge),
but can also be considerably slower (e.g. multi-processor AMD Interlagos machines). However,
there is a more pronounced benefit of multi-level parallelization in highly parallel runs.

Hybrid/heterogeneous acceleration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|Gromacs| introduces hybrid acceleration by making use of GPUs to accelerate non-bonded force
calculation. Along the :doc:`Verlet cut-off scheme <cutoff-schemes>` new non-bonded algorithms
have been developed with the aim of efficient acceleration both on CPUs and GPUs.
 
To efficiently use all compute resource available, CPU and GPU computation is done simultaneously.
Overlapping with the OpenMP multithreaded bonded force and PME long-range electrostatic calculations
on the CPU, non-bonded forces are calculated on the GPU. Multiple GPUs, both in a single node as
well as across multiple nodes, are supported using domain-decomposition. A single GPU is assigned
to the non-bonded workload of a domain, therefore, the number GPUs used has to match the number
of of MPI processes (or thread-MPI threads) the simulation is started with. That the available
CPU cores are partitioned among the processes (or thread-MPI threads) and a set of cores
with a GPU do the calculations on the respective domain.
 
With PME electrostatics, :ref:`mdrun <gmx mdrun>` supports automated CPU-GPU load-balancing by
shifting workload from the PME mesh calculations, done on the CPU, to the particle-particle
non-bonded calculations, done on the GPU. At startup a few iterations of tuning are executed
during the first 100 to 1000 MD steps. These iterations involve scaling the electrostatics cut-off
and PME grid spacing to determine the value that gives optimal CPU-GPU load balance. The cut-off
value provided using the :mdp:`rcoulomb` ``=rvdw`` :ref:`mdp` option represents the minimum
electrostatics cut-off the tuning starts with and therefore should be chosen as small as
possible (but still reasonable for the physics simulated). The Lennard-Jones cut-off ``rvdw``
is kept fixed. We don't allow scaling to shorter cut-off as we don't want to change ``rvdw``
and there would be no performance gain in the Verlet cut-off scheme.
 
While the automated CPU-GPU load balancing always attempts to find the optimal cut-off setting,
it might not always be possible to balance CPU and GPU workload. This happens when the CPU threads
finish calculating the bonded forces and PME faster than the GPU the non-bonded force calculation,
even with the shortest possible cut-off. In such cases the CPU will wait for the GPU and this
time will show up as ``Wait GPU local`` in the cycle and timing summary table at the end of the log file as shown below.

.. table::

     +--------------------------------------------------------------------------+
     |  R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G             |
     +====================+=======+=======+========+=========+==========+=======+
     | Computing:         | Nodes | Th.   | Count  | Seconds | G-Cycles | %     |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Neighbor search    | 1     | 4     | 26     | 0.145   | 1.866    | 5.2   |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Launch GPU ops.    | 1     | 4     | 501    | 0.035   | 0.448    | 1.2   |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Force              | 1     | 4     | 501    | 0.338   | 4.349    | 12.0  |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | PME mesh           | 1     | 4     | 501    | 1.365   | 17.547   | 48.5  |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Wait GPU local     | 1     | 4     | 501    | 0.162   | 2.083    | 5.8   |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | NB X/F buffer ops. | 1     | 4     | 1002   | 0.128   | 1.645    | 4.6   |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Write traj.        | 1     | 4     | 1      | 0.180   | 2.309    | 6.4   |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Update             | 1     | 4     | 501    | 0.072   | 0.924    | 2.6   |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Constraints        | 1     | 4     | 501    | 0.322   | 4.147    | 11.5  |
     +--------------------+-------+-------+--------+---------+----------+-------+
     | Rest               | 1     |                | 0.065   | 0.833    | 2.3   |
     +--------------------+-------+----------------+---------+----------+-------+
     | Total              | 1     |                | 2.811   | 36.152   | 100.0 |
     +--------------------+-------+----------------+---------+----------+-------+

 
Performance issues with hybrid acceleration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With hybrid acceleration there are different resources working on different tasks or parts of the
work. When the work load is not balanced, some resources will be idling. An extreme example is
GPU-only code such as OpenMM, where the CPU, which is always present, idles all the time.
In Gromacs we are lucky that the bonded+PME calculation work load on the CPU often roughly
matches the non-bonded work load on the GPU. But how good the balance is will depend on
your hardware and simulation setup. There are two extreme cases of imbalance:

-    Reaction-field simulations, especially with little bonded interaction, e.g. pure water.
     Here the CPU has almost nothing to do while the GPU calculates the non-bonded forces. 
     If you use multiple GPUs, you could be lucky that the hybrid non-bonded scheme, turned on by
     :ref:`mdrun <gmx mdrun>` ``-nb gpu_cpu``, is faster. In the future we plan to balance the
     non-bonded work load between GPU and CPU.

-    Parallel simulations of a solvated macro-molecule with PME. When running on many GPUs,
     the domains corresponding to the protein will have a much higher work load, as with
     GPU acceleration the bonded forces start taking a significant amount of time.
     This leads to load imbalance and performance loss. Currently there is not much
     to do about this, except for placing your molecule and choosing the domain decomposition
     such that the molecule gets divided over multiple domains. We are working on a better solution for this issue.

Separate PME nodes
^^^^^^^^^^^^^^^^^^

By default, particle-particle (PP) and PME calculations are done in the same process one after
another. As PME requires heavy global communication, this is most of the time the limiting
factor to scaling on a large number of cores. By designating a subset of nodes for PME
calculations only, performance of parallel runs can be greatly improved.

OpenMP mutithreading in PME nodes is also possible and is supported with both group and
Verlet cut-off schemes. Using multi-threading in PME can can improve performance at high
parallelization. The reason for this is that with N>1 threads the number of processes
communicating, and therefore the number of messages, is reduced by a factor of N.
But note that modern communication networks can process several messages simultaneously,
such that it could be advantages to have more processes communicating.
 
Separate PME nodes are not used at low parallelization, the switch at higher parallelization
happens automatically (at > 16 processes). The number of PME nodes is estimated by mdrun.
If the PME load is higher than the PP load, mdrun will automatically balance the load, but
this leads to additional (non-bonded) calculations. This avoids the idling of a large fraction
of the nodes; usually 3/4 of the nodes are PP nodes. But to ensure the best absolute performance
of highly parallel runs, it is advisable to tweak this number which is automated by the g_tune_pme tool.
 
The number of PME nodes can be set manually on the :ref:`mdrun <gmx mdrun>` command line using the ``-npme``
option, the number of PME threads can be specified on the command line with ``-ntomp_pme`` or
alternatively using the ``GMX_PME_NUM_THREADS`` environment variable. The latter is especially
useful when running on compute nodes with different number of cores as it enables
setting different number of PME threads on different nodes.
 
Running simulations
-------------------

Simple examples to run |Gromacs| on :ref:`single <gmx-mdrun-single-node>` or 
:ref:`multiple <gmx-mdrun-multiple-nodes>` nodes can be found on a different page.

We assume default mdrun options wherever the explicit values are not specified. Additionally, in the examples
:ref:`mdrun_mpi <gmx mdrun>` indicates a binary compiled with real MPI, and :ref:`mdrun <gmx mdrun>` describes the (default) compiled
with |Gromacs| built-in Thread-MPI. Note that all features available with MPI are also supported
with thread-MPI so whenever ``process`` or ``MPI process`` is used, these are equivalent.

Following are more advanced examples for getting optimal performance with |Gromacs| and different
parallelisation schemes.

Approaching the scaling limit
-----------------------------

There are several aspects of running a |Gromacs| simulation that are important as the number
of atoms per core approaches the current scaling limit of ~100 atoms/core.

One of these is that the use of ``constraints = all-bonds``  with P-LINCS
sets an artificial minimum on the size of domains. You should reconsider the use
of constraints to all bonds (and bear in mind possible consequences on the safe maximum for dt),
or change lincs_order and lincs_iter suitably.
