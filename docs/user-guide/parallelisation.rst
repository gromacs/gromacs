.. _gmx-parallel:

Acceleration and parallelization
================================

Here we give an overview on the parallelization and acceleration schemes employed by |Gromacs|
starting from version 4.6, which introduced several new acceleration and parallelization
features. Up to |Gromacs| 4.5, the situation was much simpler. The aim is on the one hand,
to provide an understanding of the underlying mechanism that make |Gromacs| one of the
fastest molecular dynamics packages. On the other hand, the information presented
should help choosing the appropriate parallelization options, run configuration,
as well as acceleration options in order to achieve optimal simulation performance.

Terms and definitions
---------------------

* Node

  It general refers to a single computer, either a workstation or a machine in
  a computer network. In |Gromacs| terminology a node can also refer to a process
  or thread-MPI thread in charge of a certain part of the simulation domain.

* Core, physical core, virtual core, "Bulldozer" module

  A core is the computational unit of a multi-core processor traditionally equivalent
  to a physical core. One physical core can support multiple logical cores or hardware
  threads. Modern Intel CPUs with `Hyper-threading (HT) <http://en.wikipedia.org/wiki/HyperThreading>`__
  capable of `simultaneous multithreading <http://en.wikipedia.org/wiki/Simultaneous_multi-threading>`__
  exposed to the operating system through virtual (logical) cores.
  At the same time, the `AMD Bulldozer <http://en.wikipedia.org/wiki/Amd_bulldozer>`__
  microarchitecture uses clustered multithreading
  in form of modules consisting cores quite different from the traditional physical cores.

* Accelerator, GPU
  
  `Graphics processing units (GPUs) <http://en.wikipedia.org/wiki/Graphics_processing_unit>`__
  are powerful compute-accelerators with strong
  floating point capabilities. |Gromacs| makes use of GPUs with the native GPU
  acceleration support in v4.6. The OpenMM-based acceleration, introduced in version 4.5,
  which runs entirely on GPU has been moved to contrib and is not actively supported.

* Thread-MPI, OpenMP
  Used in parallelization within a node, multithreading enables efficient use of multicore
  CPUs. Multithreading was first introduced in |Gromacs| 4.5 based on thread-MPI library
  which provides a threading-based MPI implementation. OpenMP-based multithreading is
  supported with |Gromacs| 4.6 and can be combined with (thread-)MPI parallelization.

* Accelerated code, SSE, AVX, CUDA
  
  To achieve high computational efficiency, |Gromacs| uses both CPU- and GPU-based
  acceleration. The most compute-intensive parts of the code are implemented as
  accelerated compute kernels for CPU using SSE or AVX and for GPUs using CUDA. 

* Heterogeneous parallelization (CPU + accelerator)
  
  A hybrid or heterogeneous parallelization makes use of multiple different computational
  units, typically CPUs and GPUs. With the native GPU acceleration support, |Gromacs| 4.6
  introduces hybrid parallelization.

* Hybrid/multi-level parallelization

  Consists of the use of multiple parallelization schemes on different hardware levels
  typically separating intra- and inter-node parallelization. |Gromacs| uses OpenMP
  multithreading for intra-node multicore-targeting parallelization and MPI for inter-node.


Acceleration
------------

CPU acceleration: SSE, AVX, etc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before Gromacs 4.6, there were mainly assembly non-bonded kernels and some optimized C and
Fortran non-bonded kernels, as well as some SSE intrinsics in PME and the implicit
solvation functionality. In Gromacs 4.6 the assembly and Fortran non-bonded kernels have been
removed. These have been replaced by three levels of non-bonded kernels: reference or generic
kernels, optimized ``plain-C`` kernels and SIMD intrinsic accelerated kernels. Other compute
intensive parts of the code, mainly PME, bonded force calculation, and neighbour searching
also employ SIMD intrinsic acceleration.

Unlike the old assembly kernels, the new SIMD intrinsic code is compiled by the compiler.
Technically, it is possible to compile different levels of acceleration into one binary,
but this is difficult to manage with acceleration in many parts of the code. Thus, you need
to configure and compile |Gromacs| with a single target hardware acceleration which corresponds
to a SIMD instruction set. By default, the build system will detect the highest supported
acceleration of the host where the compilation is carried out. For cross-compiling for 
a machine with different highest SIMD instructions set, in order to set the target acceleration,
the ``-DGMX_CPU_ACCELERATION`` CMake option can be used. For best performance always pick the highest
(latest) SIMD instruction set supported by the target architecture (and |Gromacs|). To use a single
installation on multiple different machines, it is convenient to compile the analysis tools with
the lowest common SIMD instruction set (as these rely little on SIMD acceleration), but for best
performance mdrun should be compiled separately for each machine.

Currently the supported acceleration options are: none, SSE2, SSE4.1, AVX-128-FMA
(AMD Bulldozer + Piledriver), AVX-256 (Intel Sandy+Ivy Bridge) and AVX2 (Intel Haswell/Haswell-E,Skylake).
We will add Blue Gene P and/or Q. On x86, the performance difference between SSE2 and SSE4.1 is minor.
All other, higher acceleration  differences are significant. Another effect of switching to intrinsics
is that the choice of compiler now affects the performance. On x86 we advice the GNU compilers (gcc) version
4.7 or later or Intel Compilers version 12 or later. Different parts of the code on different CPUs can
see performance differences of up to 10% between these two compilers, in either direction. At the time
of writing, in most of our benchmarks we observed gcc 4.7/4.8 to generate faster code.

GPU acceleration
^^^^^^^^^^^^^^^^

|Gromacs| 4.5 introduced the first version of GPU acceleration based on the
`OpenMM library <https://simtk.org/home/openmm>`. This version executed the entire simulation on the
GPU and doesn't use the CPU resources for anything but input-ouput. While this approach avoids
the CPU-GPU communication bottleneck, it only supports a
`small subset <http://www.gromacs.org/Documentation/Installation_Instructions_4.5/GROMACS-OpenMM#Supported_features>`
of all |Gromacs| features and delivers substantial speedup compared to CPU runs only in case of
`implicit solvent simulations <http://www.gromacs.org/Documentation/Installation_Instructions_4.5/GROMACS-OpenMM#Benchmark_results.3a_GROMACS_CPU_vs_GPU>`.
 
With |Gromacs| 4.6, native :ref:`GPU acceleration <gmx-mdrun-on-gpu>` support is introduced.
The most compute-intensive part of simulations, the non-bonded force calculation can be
offloaded a GPU and carried out simultaneously with CPU calculations of bonded forces and
PME eletrostatics. Native GPU acceleration is supported with the :doc:`verlet cut-off scheme <cutoff-schemes>`
(not with the group scheme) with PME, reaction-field, and plain cut-off electrostatics.    
 
The native non-bonded GPU kernels are implemented in CUDA and run on any NVIDIA hardware with
compute capability 2.0 or higher, that is GPUs with ``Fermi`` and ``Kepler`` chips or higher. Support is
not limited to high-end cards and professional cards like Tesla and Quadro, GeForce cards work equally well.
Although low-end GPUs (e.g. GeForce GTX 630) will work, typically at least a mid-class consumer GPU is
needed to achieve speedup compared to CPU-only runs on a recent processor. For optimal performance with
multiple GPUs, especially in multi-node runs, it is best to use identical hardware as balancing
the load between different GPU is not possible.
 
The native GPU acceleration can be turned on or off, either at run-time using the
:ref:`mdrun <gmx mdrun>` ``-nb`` option, or at configuration time using the ``GMX_GPU`` CMake variable.
 
From |Gromacs| 5.0+, native GPU acceleration supports now both CUDA and OpenCL. With CUDA, it is also
optimized on Maxwell architectures (CUDA Compute Capability 5.0/5.2). OpenCL currently works well only
in MacOS X and AMD GPUs. The Verlet scheme still includes only analytical non-bonded Van der Waals interactions.
Soon tabulated potentials for non-bonded generic, Coulomb and Van der Waals will be fully supported in
CUDA from |Gromacs| version 5.1.x. If you are interested in downloading our first working CUDA implementation,
please use `this patch <https://gerrit.gromacs.org/#/c/5046/>`. Currently, only the LJ supports multiple
atom types, while the generic table supports only a single table.

Parallelization schemes
-----------------------

|Gromacs|, being performance-oriented, has a strong focus on efficient parallelization. As of version
4.6, there are multiple parallelization schemes available, therefore a simulation can be run on a
given hardware with different choices of run configuration. Here we describe the different schemes
employed in |Gromacs| 4.6, highlighting the differences and providing a guide for running efficient simulations.

MPI
^^^

Parallelization based on MPI has been part of |Gromacs| from the early versions hence is compatible
with the majority of MD algorithms. At the heart of the MPI parallelization is the neutral-territory
:ref:`domain decomposition <gmx-domain-decomp>` which supports fully automatic dynamic load balancing.
Particle decomposition is also supported with MPI.

To parallelize simulations across multiple machines (e.g. nodes of a cluster) 
:ref:`mdrun <gmx mdrun>` needs to be compiled with MPI which can be enbled using the ``GMX_MPI`` CMake variable.
 
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
``-nt`` and ``-ntmpi`` options; in 4.5 only the former is supported as thread-MPI is the
only means of multi-threading, but in 4.6 ``-nt`` represents the total number of threads
to be used (which can be a mix of thread-MPI and OpenMP threads with the :doc:`verlet scheme <cutoff-schemes>`).
Note that in version 4.5 and 4.6, if the number of threads :ref:`mdrun <gmx mdrun>`
uses is equal with the total number of cores, each thread gets locked to ``its`` core.
 
 
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
 
With |Gromacs| 4.6, OpenMP multithreading is supported in :ref:`mdrun <gmx mdrun>`
and combined with MPI (or thread-MPI) it enables multi-level and heterogeneous parallelization.
With the :doc:`verlet cut-off scheme <cutoff-schemes>` full OpenMP multithreading support is implemented,
but the group scheme currently only supports OpenMP threading for PME. 
 
OpenMP is enabled by default in |Gromacs| 4.6 and can be turned on/off at configure time with
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

|Gromacs| 4.6 introduces hybrid acceleration by making use of GPUs to accelerate non-bonded force
calculation. Along the :doc:`verlet cut-off scheme <cutoff-schemes>` new non-bonded algorithms
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
and there would be no performance gain in the verlet cut-off scheme.
 
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

Using separate PME nodes has been possible since |Gromacs| 4.0. With version 4.6
OpenMP mutithreading in PME nodes is also possible and is supported with both group and
verlet cut-off schemes. Using multi-threading in PME can can improve performance at high
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

Pinning threads to physical cores
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, thread-MPI and OpenMP parallelization fill up all cores in the machine. When all
cores are used, mdrun will pin the threads to specific cores (also known as setting the thread
affinities for the cores), unless it detects this has already been done (e.g. by MPI or OpenMP).
This stops the operating system kernel from moving |Gromacs| processes between cores, which it might
otherwise have done in response to non-\ |Gromacs| processes being run on the machine. Being
able to move a |Gromacs| process when all cores have |Gromacs| processes is generally more
wasteful than waiting for the old core to become free.

If you want optimal performance when not using all cores, you need to use :ref:`mdrun <gmx mdrun>` ``-pin on``.
This is particularly true if your hardware is heterogeneous and not evenly divisible
(e.g. 3 GPUs on a node with four quad-core sockets).

If you want to run multiple jobs on the same compute node, you need to limit the number of cores
used and if you want good performance, pin different jobs to different cores. The :ref:`mdrun <gmx mdrun>`
option ``-nt`` sets the total number of threads for an :ref:`mdrun <gmx mdrun>` job. The ``-pinoffset``
option sets a pinning offset, counted in logical numbers of cores. For example,
running 2 :ref:`mdrun <gmx mdrun>` jobs on an Intel CPU with 6 physical cores with hyper-threading
g (supporting 2 threads) can be achieved with:

::

    mdrun -nt 6 -pin on -pinoffset 0
    mdrun -nt 6 -pin on -pinoffset 3

Multi-level parallelization: MPI/thread-MPI + OpenMP

Combining MPI/thread-MPI with OpenMP has a considerable overhead. Therefore, at the moment, the
multi-level parallelization will surpass the (thread-)MPI-only parallelization only in case of
highly parallel runs and/or with a slow network. We refer to the verlet scheme unless explicitly
stated as only this scheme has full OpenMP support. 

Launching M MPI processes with N OpenMP threads each:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::
    mdrun -ntmpi M -ntomp N
    mpirun -np M mdrun_mpi -ntomp N

But as you usually want to use all available hardware, the ``-ntomp`` option can be omitted:

::
    mdrun -ntmpi M
    mpirun -np M mdrun_mpi

Note that for good performance on multi-socket servers, groups of OpenMP threads belonging to an
MPI process/thread-MPI thread should run on the same CPU/socket. This requires that the number of processes
is a multiple of the number of CPUs/sockets in the respective machine and the number of cores per CPU
is divisible by the number of threads per process. E.g. on a dual 6-core machine N=6, M=2
or N=3, M=4 should run more efficiently than N=4 and M=3.

Running separate, multi-threaded PME nodes is supported in both cut-off schemes. To set the number
of threads for PME only independently from the number of threads in the rest of the code, there
is the ``-ntomp`` option. While with the verlet scheme it is mandatory to always set the global
number of threads (``-ntomp``) if the number of PME threads is set, with the group scheme it is enough ``-ntomp_pme``.

Examples:

::
    mpirun -np NP_tot mdrun_mpi -npme NP_pme -ntomp NT

will run NP_tot processes out of which NP_pme dedicated for PME using NT threads for both PP and PME, while 

::
    mpirun -np NP_tot mdrun_mpi -npme NP_pme -ntomp NT -ntomp_pme NT_pme

will use NT threads in PP nodes and NT_pme threads in PME nodes.

 
Heterogenous parallelization: using GPUs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using GPU acceleration is pretty much as simple as compiling mdrun with the CMake variable ``GMX_GPU=ON``
and using a tpr file with the Verlet scheme on a machine with supported GPU(s). Therefore, the above
instructions regarding OpenMP and MPI/thread-MPI + OpenMP runs apply to GPU accelerated runs too.
The only restriction with GPU runs is that the current parallelization scheme uses domain-decomposition
to utilize multiple GPUs by assigning the computation of non-bonded forces in a domain to a GPU on the same
physical node. Therefore, the number of GPUs used determines the domain-decomposition required, e.g with
four GPUs at least four-way decomposition is needed with four particle-particle ranks. Hence, the number of
ranks mdrun is started with (be it PP+PME or PP-only with separate PME) has to be equal with (or multiple of)
the number of available GPUs. Consequently, you need to make sure to start a number of MPI ranks
that is a multiple of the number of GPUs intended to be used. With thread-MPI the number of MPI threads
is automatically set to the number of compatible GPUs (note that this could include slow GPUs).

For instance, with an 8-core machine with two GPUs the launch command with thread-MPI can be as simple as:

::
    mdrun

as in this case we detect the two GPUs, start two MPI threads with one GPU assigned to each. This is equivalent with the following:

::
    mdrun -ntmpi 2
    mdrun -ntmpi 2 -ntomp 4 #2 x 4 threads = 8 threads

and with MPI:

::
    mpirun -np 2 mdrun_mpi

GPUs are assigned to PP ranks within the same physical node in a sequential order, that is GPU 0 to
the (thread-)MPI rank 0, GPU 1 to rank 1. In order to manually specify which GPU(s) to be used by
:ref:`mdrun <gmx mdrun>`, the respective device ID(s) can be passed with the ``-gpu_id XYZ``
command line option or with the ``GMX_GPU_ID=XYZ`` environment variable. Here, XYZ is a
sequence of digits representing the numeric ID-s of available GPUs (the numbering starts from 0).
The environment variable is particularly useful when running on multiple compute nodes with different GPU configurations.

Taking the above example of 8-core machine with two compatible GPUs, we can manually specify the GPUs
and get the same launch configuration as in the above examples by:

::
    mdrun -ntmpi 2 -ntomp 4 -gpu_id 01

Now, let's assume that the first GPU in our 8-core system is a slow one used only for driving the
display. In this case we typically want to avoid using this GPU for computation which can be achieved by running:

::
    mdrun -gpu_id 1 # the first device, GPU0, will not be used

Note that in this example  mdrun will know that we intend to use only a single GPU which requires a
single domain (i.e no domain-decomposition) and therefore will start a single MPI rank with OpenMP 8 threads. 
If we add a third GPU for compute use we have to modify the above command to:

::
    mdrun -gpu_id 12 # skip GPU0, use the 2nd and 3rd device

If we run across 5 nodes of a cluster, with one PP rank per node and each node with one GPU per node, we could use

::
    mpirun -np 5 mdrun_mpi -gpu_id 0 # Use GPU zero; in this case, specifying -gpu_id is optional

If we wanted two PP ranks per node on a 5-node machine, we could use

::
    mpirun -np 10 mdrun_mpi -gpu_id 01

Currently, the automation of GPU to process assignment is fairly simplistic, GPUs will be automatically
assigned sequentially to threads/processes meaning that the (PP/PP+PME) process IDs within a machine will
match the GPU ID. Although this scheme works well in the majority of cases, it does not take into account
locality (on the PCI-E bus) and the performance of each GPU, each GPU will be assumed to have the same
performance. Additionally, when more GPUs are available than processes/threads started (by specifying
``mpirun -np N`` or ``-ntmpi N``), :ref:`mdrun <gmx mdrun>` does a ``naive`` choice and will use the
first N GPUs rather than the fastest ones. In this case the fast GPUs intended to be used need to be
manually specified (using ``-gpu_id`` or ``GMX_GPU_ID``) skipping the GPU devices not intended to be used.

 
Multiple MPI ranks per GPU
^^^^^^^^^^^^^^^^^^^^^^^^^^

As explained earlier, when using GPU acceleration, the short-range non-bonded forces are calculated on the
GPU while the CPU calculated bonded forces and Ewald long-range electrostatics (with PME). CPU cores working
in parallel with the GPU need to belong to the same ``team`` of OpenMP threads, hence to the same MPI rank.
Therefore, the number of GPUs in a compute node will typically determine the number of (PP) MPI ranks
needed, hence the number of threads per rank. However, the OpenMP multi-threaded parallelization is rather
sensitive and it often does not scale well to large number of threads, especially with teams of threads in
a ranks spanning across CPUs/NUMA regions. The potential slowdowns get more pronounced when running in
parallel on multiple compute nodes. In these cases, to address the bottleneck caused by multi-threading
inefficiencies, it can be advantageous to reduce the number of OpenMP threads per rank. However, to not
leave cores empty, this requires using more MPI ranks, hence more PP ranks, and therefore ranks will have
to share GPUs. GPU sharing is possible by passing a GPU ID to :ref:`mdrun <gmx mdrun>` multiple times, e.g
``-gpu_id 0011`` will allow the first two PP ranks in a compute node to use GPU0 and the third and fourth GPU1.

For instance, given a dual-socket AMD Opteron machine with two 6-core CPUs and a fast GPU, like a
GeForce GTX680 or Tesla K20, simply starting mdrun with the default launch configuration will
lead to a run equivalent with the following:

::
    mpirun -np 1 mdrun -ntomp 12 -gpu_id 0 # equivalent with the default launch config mdrun will use

This means that a single MPI rank with 12 OpenMP threads will be used together with the GPU. Such a
configuration that runs many OpenMP threads per MPI rank will often be hampered by inefficient multithreading,
e.g. on AMD Opteron processors with two NUMA regions as threads will communicate through bus linking the
two dies on a chip. To address this, we can try to run multiple MPI ranks per GPU with fewer threads each,
e.g. two ranks with 6 threads or four ranks with 3 threads each:

::
    mpirun -np 2 mdrun -ntomp 6 -gpu_id 00 # two ranks sharing GPU0
    mpirun -np 4 mdrun -ntomp 3 -gpu_id 0000 # four ranks sharing GPU0

Of course it is possible to do the same using multiple compute nodes, e.g. on 5 identical nodes of a cluter: 

::
    mpirun -np 10 mdrun -ntomp 6 -gpu_id 00 # two ranks sharing GPU0
    mpirun -np 20 mdrun -ntomp 3 -gpu_id 0000 # four ranks sharing GPU0

Remarks:

-    On Intel machines, especially if running only a single compute-node, as the OpenMP
     multi-threading bottlenecks are less severe than on AMD, it can be faster to not use
     domain-decomposition (which itself imposes a certain overhead), but instead run OpenMP
     threads across CPUs (like in the first example). However, on newer clusters with Sandy Bridge
     or Ivy Bridge processors with 10-12 cores it is most of the time more advantageous
     to also share a GPU among multiple PP ranks.

-    In versions 4.6-4.6.4, the measured and reported domain-decomposition load imbalance was
     usually incorrect when sharing GPUs, and tuning off load balancing (``-dlb no``)
     could actually improve performance in some cases. Fixed in 4.6.5.
     

Using multi-simulations and GPUs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using :ref:`mdrun <gmx mdrun>` ``-multi`` to run multiple simulations in one call of mpirun
(e.g. for REMD) is also supported with GPUs. There still needs to be a mapping of PP MPI
ranks to GPU ids, but those PP ranks do not all have to come from the same component simulation.
The mapping of MPI ranks into component simulations is distinct from the mapping of PP MPI ranks
to GPUs. There are degenerate cases where you will not need to specify ``-gpu_id``. For
example, on a machine with 4 physical nodes, with 2 GPUs per physical node with MPI configured
to produce 16 MPI processes per physical node, you can use

::
    pirun -np 64 mdrun-mpi-gpu -multi 4 -gpu_id 0000000011111111

to run a four-component multi-simulation.

Note that it is most often advantageous to run multiple independent simulations (either part of a
multi-sim or not) on a single GPU. In the single simulation per GPU case, the GPU utilization is
limited to the amount of possible overlap between CPU and GPU computation during a time-step.
In contrast, independent simulations do not need to synchronize every time-step and can significantly
increase the overall GPU utilization. As a consequence, multiple independent runs (part of a
multi-sim or not) using the same GPU will most often lead to considerably higher aggregate
simulation speed when run simultaneously, compared to running them in a sequence.

Approaching the scaling limit
-----------------------------

There are several aspects of running a |Gromacs| simulation that are important as the number
of atoms per core approaches the current scaling limit of ~100 atoms/core.

One of these is that the use of

::
    constraints = all-bonds

with P-LINCS sets an artificial minimum on the size of domains. You should reconsider the use
of constraints to all bonds (and bear in mind possible consequences on the safe maximum for dt),
or change lincs_order and lincs_iter suitably.
