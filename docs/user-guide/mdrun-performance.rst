.. _gmx-performance:

Getting good performance from :ref:`mdrun <gmx mdrun>`
======================================================

Here we give an overview on the parallelization and acceleration schemes employed by |Gromacs|.
The aim is to provide an understanding of the underlying mechanisms that make |Gromacs| one of the
fastest molecular dynamics packages. The information presented
should help choosing appropriate parallelization options, run configuration,
as well as acceleration options to achieve optimal simulation performance.


The |Gromacs| build system and the :ref:`gmx mdrun` tool have a lot of built-in
and configurable intelligence to detect your hardware and make pretty
effective use of it. For a lot of casual and serious use of
:ref:`gmx mdrun`, the automatic machinery works well enough. But to get the
most from your hardware to maximize your scientific quality, read on!

Hardware background information
-------------------------------
Modern computer hardware is complex and heterogeneous, so we need to
discuss a little bit of background information and set up some
definitions. Experienced HPC users can skip this section.

.. glossary::

    core
        A hardware compute unit that actually executes
        instructions. There is normally more than one core in a
        processor, often many more.

    cache
        A special kind of memory local to core(s) that is much faster
        to access than main memory, kind of like the top of a human's
        desk, compared to their filing cabinet. There are often
        several layers of caches associated with a core.

    socket
        A group of cores that share some kind of locality, such as a
        shared cache. This makes it more efficient to spread
        computational work over cores within a socket than over cores
        in different sockets. Modern processors often have more than
        one socket.

    node
        A group of sockets that share coarser-level locality, such as
        shared access to the same memory without requiring any network
        hardware. A normal laptop or desktop computer is a node. A
        node is often the smallest amount of a large compute cluster
        that a user can request to use.

    thread
        A stream of instructions for a core to execute. There are many
        different programming abstractions that create and manage
        spreading computation over multiple threads, such as OpenMP,
        pthreads, winthreads, CUDA, OpenCL, and OpenACC. Some kinds of
        hardware can map more than one software thread to a core; on
        Intel x86 processors this is called "hyper-threading", while
        the more general concept is often called SMT for
        "simultaneous multi-threading". IBM Power8 can for instance use
        up to 8 hardware threads per core.
        This feature can usually be enabled or disabled either in
        the hardware BIOS or through a setting in the Linux operating
        system. |Gromacs| can typically make use of this, for a moderate
        free performance boost. In most cases it will be
        enabled by default e.g. on new x86 processors, but in some cases
        the system administrators might have disabled it. If that is the
        case, ask if they can re-enable it for you. If you are not sure
        if it is enabled, check the output of the CPU information in
        the log file and compare with CPU specifications you find online.

    thread affinity (pinning)
        By default, most operating systems allow software threads to migrate
        between cores (or hardware threads) to help automatically balance
        workload. However, the performance of :ref:`gmx mdrun` can deteriorate
        if this is permitted and will degrade dramatically especially when
        relying on multi-threading within a rank. To avoid this,
        :ref:`gmx mdrun` will by default
        set the affinity of its threads to individual cores/hardware threads,
        unless the user or software environment has already done so
        (or not the entire node is used for the run, i.e. there is potential
        for node sharing).
        Setting thread affinity is sometimes called thread "pinning".

    MPI (Message Passing Interface)
        The dominant multi-node parallelization-scheme, which provides
        a standardized language in which programs can be written that
        work across more than one node.

    rank
        In MPI, a rank is the smallest grouping of hardware used in
        the multi-node parallelization scheme. That grouping can be
        controlled by the user, and might correspond to a core, a
        socket, a node, or a group of nodes. The best choice varies
        with the hardware, software and compute task. Sometimes an MPI
        rank is called an MPI process.

    GPU
        A graphics processing unit, which is often faster and more
        efficient than conventional processors for particular kinds of
        compute workloads. A GPU is always associated with a
        particular node, and often a particular socket within that
        node.

    OpenMP
        A standardized technique supported by many compilers to share
        a compute workload over multiple cores. Often combined with
        MPI to achieve hybrid MPI/OpenMP parallelism.

    CUDA
        A proprietary parallel computing framework and API developed by NVIDIA
        that allows targeting their accelerator hardware.
        |Gromacs| uses CUDA for GPU acceleration support with NVIDIA hardware.

    OpenCL
        An open standard-based parallel computing framework that consists
        of a C99-based compiler and a programming API for targeting heterogeneous
        and accelerator hardware. |Gromacs| uses OpenCL for GPU acceleration
        on AMD devices (both GPUs and APUs), Intel integrated GPUs, and Apple
        Silicon integrated GPUs; some NVIDIA hardware is also supported.
        In |Gromacs|, OpenCL has been deprecated in favor of SYCL.

    SYCL
        An open standard based on C++17 for targeting heterogeneous systems.
        SYCL has several implementations, of which |Gromacs| supports two:
        `Intel oneAPI DPC++`_ and AdaptiveCpp_. |Gromacs| uses SYCL for GPU acceleration
        on AMD and Intel GPUs. There is experimental support for NVIDIA GPUs too.

    SIMD
        A type of CPU instruction by which modern CPU cores can execute multiple
        floating-point instructions in a single cycle.


Work distribution by parallelization in |Gromacs|
-------------------------------------------------

The algorithms in :ref:`gmx mdrun` and their implementations are most relevant
when choosing how to make good use of the hardware. For details,
see the :ref:`Reference Manual <gmx-reference-manual-rst>`. The most important of these are

.. _gmx-domain-decomp:

.. glossary::

    Domain Decomposition
        The domain decomposition (DD) algorithm decomposes the
        (short-ranged) component of the non-bonded interactions into
        domains that share spatial locality, which permits the use of
        efficient algorithms. Each domain handles all of the
        particle-particle (PP) interactions for its members, and is
        mapped to a single MPI rank. Within a PP rank, OpenMP threads
        can share the workload, and some work can be offloaded to a
        GPU. The PP rank also handles any bonded interactions for the
        members of its domain. A GPU may perform work for more than
        one PP rank, but it is normally most efficient to use a single
        PP rank per GPU and for that rank to have thousands of
        particles. When the work of a PP rank is done on the CPU,
        :ref:`mdrun <gmx mdrun>` will make extensive use of the SIMD
        capabilities of the core. There are various
        :ref:`command-line options <controlling-the-domain-decomposition-algorithm>`
        to control the behaviour of the DD algorithm.

    Particle-mesh Ewald
        The particle-mesh Ewald (PME) algorithm treats the long-ranged
        component of the non-bonded interactions (Coulomb and possibly also
        Lennard-Jones).  Either all, or just a subset of ranks may
        participate in the work for computing the long-ranged component
        (often inaccurately called simply the "PME"
        component). Because the algorithm uses a 3D FFT that requires
        global communication, its parallel efficiency gets worse as more ranks
        participate, which can mean it is fastest to use just a subset
        of ranks (e.g.  one-quarter to one-half of the ranks). If
        there are separate PME ranks, then the remaining ranks handle
        the PP work. Otherwise, all ranks do both PP and PME work.

Parallelization schemes
-----------------------

|Gromacs|, being performance-oriented, has a strong focus on efficient parallelization.
There are multiple parallelization schemes available, therefore a simulation can be run on a
given hardware with different choices of run configuration.

.. _intra-core-parallelization:

Intra-core parallelization via SIMD: SSE, AVX, etc.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One level of performance improvement available in |Gromacs| is through the use of
``Single Instruction Multiple Data (SIMD)`` instructions. In detail information
for those can be found under :ref:`SIMD support <gmx-simd-support>` in the installation
guide.

In |Gromacs|, SIMD instructions are used to parallelize the parts of the code with
the highest impact on performance (nonbonded and bonded force calculation,
PME and neighbour searching), through the use of hardware specific SIMD kernels.
Those form one of the three levels of non-bonded kernels that are available: reference or generic
kernels (slow but useful for producing reference values for testing),
optimized plain-C kernels (can be used cross-platform but still slow)
and SIMD intrinsics accelerated kernels.

The SIMD intrinsic code is compiled by the compiler.
Technically, it is possible to compile different levels of acceleration into one binary,
but this is difficult to manage with acceleration in many parts of the code.
Thus, you need to configure and compile |Gromacs| for the SIMD capabilities of the target CPU.
By default, the build system will detect the highest supported
acceleration of the host where the compilation is carried out. For cross-compiling for
a machine with a different highest SIMD instructions set, in order to set the target acceleration,
the ``-DGMX_SIMD`` CMake option can be used.
To use a single
installation on multiple different machines, it is convenient to compile the analysis tools with
the lowest common SIMD instruction set (as these rely little on SIMD acceleration), but for best
performance :ref:`mdrun <gmx mdrun>` should be compiled be compiled separately with the
highest (latest) ``native`` SIMD instruction set of the target architecture (supported by |Gromacs|).

Recent Intel CPU architectures bring tradeoffs between the maximum clock frequency of the
CPU (ie. its speed), and the width of the SIMD instructions it executes (ie its throughput
at a given speed). In particular, the Intel ``Skylake`` and ``Cascade Lake`` processors
(e.g. Xeon SP Gold/Platinum), can offer better throughput when using narrower SIMD because
of the better clock frequency available. Consider building :ref:`mdrun <gmx mdrun>`
configured with ``GMX_SIMD=AVX2_256`` instead of ``GMX_SIMD=AVX512`` for better
performance in GPU accelerated or highly parallel MPI runs.

Some of the latest ARM based CPU, such as the Fujitsu A64fx, support the Scalable Vector Extensions (SVE).
Though SVE can be used to generate fairly efficient Vector Length Agnostic (VLA) code,
this is not a good fit for |Gromacs| (as the SIMD vector length assumed to be known at
CMake time). Consequently, the SVE vector length must be fixed at CMake time. The default
is to automatically detect the default vector length at CMake time
(via the ``/proc/sys/abi/sve_default_vector_length`` pseudo-file, and this can be changed by
configuring with ``GMX_SIMD_ARM_SVE_LENGTH=<len>``.
The supported vector lengths are 128, 256, 512 and 1024. Since the SIMD short-range non-bonded kernels
only support up to 16 floating point numbers per SIMD vector, 1024 bits vector length is only
valid in double precision (e.g. ``-DGMX_DOUBLE=on``).
Note that even if :ref:`mdrun <gmx mdrun>` does check the SIMD vector length at runtime, running with a different
vector length than the one used at CMake time is undefined behavior, and :ref:`mdrun <gmx mdrun>` might crash before reaching
the check (that would abort with a user-friendly error message).

Process(-or) level parallelization via OpenMP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|Gromacs| :ref:`mdrun <gmx mdrun>` supports OpenMP multithreading for all parts
of the code. OpenMP is enabled by default and
can be turned on/off at configure time with the ``GMX_OPENMP`` CMake variable
and at run-time with the ``-ntomp`` option (or the ``OMP_NUM_THREADS`` environment variable).
The OpenMP implementation is quite efficient and scales well for up to 12-24 threads on
Intel and 6-8 threads on AMD CPUs.

Node level parallelization via GPU offloading and thread-MPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Multithreading with thread-MPI
..............................

The thread-MPI library implements a subset of the MPI 1.1 specification,
based on the system threading support. Both POSIX pthreads and Windows threads are supported,
thus providing great portability to most UNIX/Linux and Windows operating systems.
Acting as a drop-in replacement for MPI, thread-MPI enables compiling and running :ref:`mdrun <gmx mdrun>`
on a single machine (i.e. not across a network) without MPI. Additionally, it not only provides a
convenient way to use computers with multicore CPU(s), but thread-MPI does in some
cases make :ref:`mdrun <gmx mdrun>` run slightly faster than with MPI.

Thread-MPI is included in the |Gromacs| source and it is the default parallelization mode,
practically rendering the serial :ref:`mdrun <gmx mdrun>` deprecated.
Compilation with thread-MPI is controlled by the ``GMX_THREAD_MPI`` CMake variable.

Thread-MPI is compatible with most :ref:`mdrun <gmx mdrun>` features and parallelization schemes,
including OpenMP, GPUs; it is not compatible with MPI and multi-simulation runs.

By default, the thread-MPI :ref:`mdrun <gmx mdrun>` will use all available cores in the machine by starting
an appropriate number of ranks or OpenMP threads to occupy all of them. The number of
ranks can be controlled using the
``-nt`` and ``-ntmpi`` options. ``-nt`` represents the total number of threads
to be used (which can be a mix of thread-MPI and OpenMP threads).

Hybrid/heterogeneous acceleration
.................................

Hybrid acceleration means distributing compute work between available CPUs and GPUs
to improve simulation performance. New non-bonded algorithms
have been developed with the aim of efficient acceleration both on CPUs and GPUs.

The most compute-intensive parts of simulations, non-bonded force calculation, as well
as possibly the PME, bonded force calculation and update and constraints can be
offloaded to GPUs and carried out simultaneously with remaining CPU work.
Native GPU acceleration is supported for the most commonly used algorithms in
|Gromacs|.
For more information about the GPU kernels, please see the :ref:`Installation guide <gmx-gpu-support>`.

The native GPU acceleration can be turned on or off, either at run-time using the
:ref:`mdrun <gmx mdrun>` ``-nb`` option, or at configuration time using the ``GMX_GPU`` CMake variable.

To efficiently use all compute resource available, CPU and GPU computation is done simultaneously.
Overlapping with the OpenMP multithreaded bonded force and PME long-range electrostatic calculations
on the CPU, non-bonded forces are calculated on the GPU. Multiple GPUs, both in a single node as
well as across multiple nodes, are supported using domain-decomposition. A single GPU is assigned
to the non-bonded workload of a domain, therefore, the number GPUs used has to match the number
of of MPI processes (or thread-MPI threads) the simulation is started with. The available
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
and there would be no performance gain.

While the automated CPU-GPU load balancing always attempts to find the optimal cut-off setting,
it might not always be possible to balance CPU and GPU workload. This happens when the CPU threads
finish calculating the bonded forces and PME faster than the GPU the non-bonded force calculation,
even with the shortest possible cut-off. In such cases the CPU will wait for the GPU and this
time will show up as ``Wait GPU NB local`` in the cycle and timing summary table at the end
of the log file.

Parallelization over multiple nodes via MPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the heart of the MPI parallelization in |Gromacs| is the neutral-territory
:ref:`domain decomposition <gmx-domain-decomp>` with dynamic load balancing.
To parallelize simulations across multiple machines (e.g. nodes of a cluster)
:ref:`mdrun <gmx mdrun>` needs to be compiled with MPI which can be enabled using the ``GMX_MPI`` CMake variable.

.. _controlling-the-domain-decomposition-algorithm:

Controlling the domain decomposition algorithm
..............................................

This section lists options that affect how the domain
decomposition algorithm decomposes the workload to the available
parallel hardware.

``-rdd``
    Can be used to set the required maximum distance for inter
    charge-group bonded interactions. Communication for two-body
    bonded interactions below the non-bonded cut-off distance always
    comes for free with the non-bonded communication. Particles beyond
    the non-bonded cut-off are only communicated when they have
    missing bonded interactions; this means that the extra cost is
    minor and nearly independent of the value of ``-rdd``. With dynamic
    load balancing, option ``-rdd`` also sets the lower limit for the
    domain decomposition cell sizes. By default ``-rdd`` is determined
    by :ref:`gmx mdrun` based on the initial coordinates. The chosen value will
    be a balance between interaction range and communication cost.

``-ddcheck``
    On by default. When inter charge-group bonded interactions are
    beyond the bonded cut-off distance, :ref:`gmx mdrun` terminates with an
    error message. For pair interactions and tabulated bonds that do
    not generate exclusions, this check can be turned off with the
    option ``-noddcheck``.

``-rcon``
    When constraints are present, option ``-rcon`` influences
    the cell size limit as well.
    Particles connected by NC constraints, where NC is the LINCS order
    plus 1, should not be beyond the smallest cell size. A error
    message is generated when this happens, and the user should change
    the decomposition or decrease the LINCS order and increase the
    number of LINCS iterations.  By default :ref:`gmx mdrun` estimates the
    minimum cell size required for P-LINCS in a conservative
    fashion. For high parallelization, it can be useful to set the
    distance required for P-LINCS with ``-rcon``.

``-dds``
    Sets the minimum allowed x, y and/or z scaling of the cells with
    dynamic load balancing. :ref:`gmx mdrun` will ensure that the cells can
    scale down by at least this factor. This option is used for the
    automated spatial decomposition (when not using ``-dd``) as well as
    for determining the number of grid pulses, which in turn sets the
    minimum allowed cell size. Under certain circumstances the value
    of ``-dds`` might need to be adjusted to account for high or low
    spatial inhomogeneity of the system.



Multi-level parallelization: MPI and OpenMP
...........................................

The multi-core trend in CPU development substantiates the need for multi-level parallelization.
Current multiprocessor machines can have 2-4 CPUs with a core count as high as 64. As the memory
and cache subsystem is lagging more and more behind the multicore evolution, this emphasizes
non-uniform memory access (NUMA) effects, which can become a performance bottleneck. At the same
time, all cores share a network interface. In a purely MPI-parallel scheme, all MPI processes
use the same network interface, and although MPI intra-node communication is generally efficient,
communication between nodes can become a limiting factor to parallelization. This is especially
pronounced in the case of highly parallel simulations with PME (which is very communication
intensive) and with ``''fat''`` nodes connected by a slow network. Multi-level parallelism aims
to address the NUMA and communication related issues by employing efficient
intra-node parallelism, typically multithreading.

Combining OpenMP with MPI creates an additional overhead
especially when running separate multi-threaded PME ranks. Depending on the architecture,
input system size, as well as other factors, MPI+OpenMP runs can be as fast and faster
already at small number of processes (e.g. multi-processor Intel Westmere or Sandy Bridge),
but can also be considerably slower (e.g. multi-processor AMD Interlagos machines). However,
there is a more pronounced benefit of multi-level parallelization in highly parallel runs.

Separate PME ranks
^^^^^^^^^^^^^^^^^^

On CPU ranks, particle-particle (PP) and PME calculations are done in the same process one after
another. As PME requires all-to-all global communication, this is most of the time the limiting
factor to scaling on a large number of cores. By designating a subset of ranks for PME
calculations only, performance of parallel runs can be greatly improved.

OpenMP multithreading in PME ranks is also possible.
Using multi-threading in PME can can improve performance at high
parallelization. The reason for this is that with N>1 threads the number of processes
communicating, and therefore the number of messages, is reduced by a factor of N.
But note that modern communication networks can process several messages simultaneously,
such that it could be advantageous to have more processes communicating.

Separate PME ranks are not used at low parallelization, the switch at higher parallelization
happens automatically (at > 16 processes). The number of PME ranks is estimated by mdrun.
If the PME load is higher than the PP load, mdrun will automatically balance the load, but
this leads to additional (non-bonded) calculations. This avoids the idling of a large fraction
of the ranks; usually 3/4 of the ranks are PP ranks. But to ensure the best absolute performance
of highly parallel runs, it is advisable to tweak this number which is automated by
the :ref:`tune_pme <gmx tune_pme>` tool.

The number of PME ranks can be set manually on the :ref:`mdrun <gmx mdrun>` command line using the ``-npme``
option, the number of PME threads can be specified on the command line with ``-ntomp_pme`` or
alternatively using the ``GMX_PME_NUM_THREADS`` environment variable. The latter is especially
useful when running on compute nodes with different number of cores as it enables
setting different number of PME threads on different nodes.

Running :ref:`mdrun <gmx mdrun>` within a single node
-----------------------------------------------------

:ref:`gmx mdrun` can be configured and compiled in several different ways that
are efficient to use within a single :term:`node`. The default configuration
using a suitable compiler will deploy a multi-level hybrid parallelism
that uses CUDA, OpenMP and the threading platform native to the
hardware. For programming convenience, in |Gromacs|, those native
threads are used to implement on a single node the same MPI scheme as
would be used between nodes, but much more efficient; this is called
thread-MPI. From a user's perspective, real MPI and thread-MPI look
almost the same, and |Gromacs| refers to MPI ranks to mean either kind,
except where noted. A real external MPI can be used for :ref:`gmx mdrun` within
a single node, but runs more slowly than the thread-MPI version.

By default, :ref:`gmx mdrun` will inspect the hardware available at run time
and do its best to make fairly efficient use of the whole node. The
log file, stdout and stderr are used to print diagnostics that
inform the user about the choices made and possible consequences.

A number of command-line parameters are available to modify the default
behavior.

``-nt``
    The total number of threads to use. The default, 0, will start as
    many threads as available cores. Whether the threads are
    thread-MPI ranks, and/or OpenMP threads within such ranks depends on
    other settings.

``-ntmpi``
    The total number of thread-MPI ranks to use. The default, 0,
    will start one rank per GPU (if present), and otherwise one rank
    per core.

``-ntomp``
    The total number of OpenMP threads per rank to start. The
    default, 0, will start one thread on each available core.
    Alternatively, :ref:`mdrun <gmx mdrun>` will honor the appropriate system
    environment variable (e.g. ``OMP_NUM_THREADS``) if set.
    Note that the maximum number of OpenMP threads (per rank) is,
    for efficiency reasons, limited to 64. While it is rarely beneficial to use
    a number of threads higher than this, the GMX_OPENMP_MAX_THREADS CMake variable
    can be used to increase the limit.

``-npme``
    The total number of ranks to dedicate to the long-ranged
    component of PME, if used. The default, -1, will dedicate ranks
    only if the total number of threads is at least 12, and will use
    around a quarter of the ranks for the long-ranged component.

``-ntomp_pme``
    When using PME with separate PME ranks,
    the total number of OpenMP threads per separate PME rank.
    The default, 0, copies the value from ``-ntomp``.

``-pin``
    Can be set to "auto," "on" or "off" to control whether
    :ref:`mdrun <gmx mdrun>` will attempt to set the affinity of threads to cores.
    Defaults to "auto," which means that if :ref:`mdrun <gmx mdrun>` detects that all the
    cores on the node are being used for :ref:`mdrun <gmx mdrun>`, then it should behave
    like "on," and attempt to set the affinities (unless they are
    already set by something else).

``-pinoffset``
    If ``-pin on``, specifies the logical core number to
    which :ref:`mdrun <gmx mdrun>` should pin the first thread. When running more than
    one instance of :ref:`mdrun <gmx mdrun>` on a node, use this option to to avoid
    pinning threads from different :ref:`mdrun <gmx mdrun>` instances to the same core.

``-pinstride``
    If ``-pin on``, specifies the stride in logical core
    numbers for the cores to which :ref:`mdrun <gmx mdrun>` should pin its threads. When
    running more than one instance of :ref:`mdrun <gmx mdrun>` on a node, use this option
    to avoid pinning threads from different :ref:`mdrun <gmx mdrun>` instances to the
    same core.  Use the default, 0, to minimize the number of threads
    per physical core - this lets :ref:`mdrun <gmx mdrun>` manage the hardware-, OS- and
    configuration-specific details of how to map logical cores to
    physical cores.

``-ddorder``
    Can be set to "interleave," "pp_pme" or "cartesian."
    Defaults to "interleave," which means that any separate PME ranks
    will be mapped to MPI ranks in an order like PP, PP, PME, PP, PP,
    PME, etc. This generally makes the best use of the available
    hardware. "pp_pme" maps all PP ranks first, then all PME
    ranks. "cartesian" is a special-purpose mapping generally useful
    only on special torus networks with accelerated global
    communication for Cartesian communicators. Has no effect if there
    are no separate PME ranks.

``-nb``
    Used to set where to execute the short-range non-bonded interactions.
    Can be set to "auto", "cpu", "gpu."
    Defaults to "auto," which uses a compatible GPU if available.
    Setting "cpu" requires that no GPU is used. Setting "gpu" requires
    that a compatible GPU is available and will be used.

``-pme``
    Used to set where to execute the long-range non-bonded interactions.
    Can be set to "auto", "cpu", "gpu."
    Defaults to "auto," which uses a compatible GPU if available.
    Setting "gpu" requires that a compatible GPU is available.
    Multiple PME ranks are not supported with PME on GPU, so if a GPU is used
    for the PME calculation -npme must be set to 1.

``-bonded``
    Used to set where to execute the bonded interactions that are part of the
    PP workload for a domain.
    Can be set to "auto", "cpu", "gpu."
    Defaults to "auto," which uses a compatible CUDA GPU only when one
    is available, a GPU is handling short-ranged interactions, and the
    CPU is handling long-ranged interaction work (electrostatic or
    LJ). The work for the bonded interactions takes place on the same
    GPU as the short-ranged interactions, and cannot be independently
    assigned.
    Setting "gpu" requires that a compatible GPU is available and will
    be used.

``-update``
    Used to set where to execute update and constraints, when present.
    Can be set to "auto", "cpu", "gpu."
    Defaults to "auto," which currently always uses the CPU.
    Setting "gpu" requires that a compatible CUDA GPU is available,
    the simulation uses a single rank.
    Update and constraints on a GPU is currently not supported
    with mass and constraints free-energy perturbation, domain
    decomposition, virtual sites, Ewald surface correction,
    replica exchange, constraint pulling, orientation restraints
    and computational electrophysiology.

``-gpu_id``
    A string that specifies the ID numbers of the GPUs that
    are available to be used by ranks on each node. For example,
    "12" specifies that the GPUs with IDs 1 and 2 (as reported
    by the GPU runtime) can be used by :ref:`mdrun <gmx mdrun>`. This is useful
    when sharing a node with other computations, or if a GPU that
    is dedicated to a display should not be used by |Gromacs|.
    Without specifying this parameter, :ref:`mdrun <gmx mdrun>`
    will utilize all GPUs. When many GPUs are
    present, a comma may be used to separate the IDs, so
    "12,13" would make GPUs 12 and 13 available to :ref:`mdrun <gmx mdrun>`.
    It could be necessary to use different GPUs on different
    nodes of a simulation, in which case the environment
    variable ``GMX_GPU_ID`` can be set differently for the ranks
    on different nodes to achieve that result.
    In |Gromacs| versions preceding 2018 this parameter used to
    specify both GPU availability and GPU task assignment.
    The latter is now done with the ``-gputasks`` parameter.

``-gputasks``
    A string that specifies the ID numbers of the GPUs to be
    used by corresponding GPU tasks on this node. For example,
    "0011" specifies that the first two GPU tasks will use GPU 0,
    and the other two use GPU 1. When using this option, the
    number of ranks must be known to :ref:`mdrun <gmx mdrun>`, as well as where
    tasks of different types should be run, such as by using
    ``-nb gpu`` - only the tasks which are set to run on GPUs
    count for parsing the mapping. See `Assigning tasks to GPUs`_
    for more details. Note that ``-gpu_id`` and
    ``-gputasks`` can not be used at the same time!
    In |Gromacs| versions preceding 2018 only a single type
    of GPU task ("PP") could be run on any rank. Now that there is some
    support for running PME on GPUs, the number of GPU tasks
    (and the number of GPU IDs expected in the ``-gputasks`` string)
    can actually be 3 for a single-rank simulation. The IDs
    still have to be the same in this case, as using multiple GPUs
    per single rank is not yet implemented.
    The order of GPU tasks per rank in the string is PP first,
    PME second. The order of ranks with different kinds of GPU tasks
    is the same by default, but can be influenced with the ``-ddorder``
    option and gets quite complex when using multiple nodes.
    Note that the bonded interactions for a PP task may
    run on the same GPU as the short-ranged work, or on the CPU,
    which can be controlled with the ``-bonded`` flag.
    The GPU task assignment (whether manually set, or automated),
    will be reported in the :ref:`mdrun <gmx mdrun>` output on
    the first physical node of the simulation. For example:

    ::

      gmx mdrun -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4

    will produce the following output in the log file/terminal:

    ::

      On host tcbl14 2 GPUs selected for this run.
      Mapping of GPU IDs to the 4 GPU tasks in the 4 ranks on this node:
      PP:0,PP:0,PP:0,PME:1

    In this case, 3 ranks are set by user to compute PP work
    on GPU 0, and 1 rank to compute PME on GPU 1.
    The detailed indexing of the GPUs is also reported in the log file.

    For more information about GPU tasks, please refer to
    :ref:`Types of GPU tasks<gmx-gpu-tasks>`.

``-pmefft``
    Allows choosing whether to execute the 3D FFT computation on a CPU or GPU.
    Can be set to "auto", "cpu", "gpu.".
    When PME is offloaded to a GPU ``-pmefft gpu`` is the default,
    and the entire PME calculation is executed on the GPU. However,
    in some cases, e.g. with a relatively slow or older generation GPU
    combined with fast CPU cores in a run, moving some work off of the GPU
    back to the CPU by computing FFTs on the CPU can improve performance.

.. _gmx-mdrun-single-node:

Examples for :ref:`mdrun <gmx mdrun>` on one node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    gmx mdrun

Starts :ref:`mdrun <gmx mdrun>` using all the available resources. :ref:`mdrun <gmx mdrun>`
will automatically choose a fairly efficient division
into thread-MPI ranks, OpenMP threads and assign work
to compatible GPUs. Details will vary with hardware
and the kind of simulation being run.

::

    gmx mdrun -nt 8

Starts :ref:`mdrun <gmx mdrun>` using 8 threads, which might be thread-MPI
or OpenMP threads depending on hardware and the kind
of simulation being run.

::

    gmx mdrun -ntmpi 2 -ntomp 4

Starts :ref:`mdrun <gmx mdrun>` using eight total threads, with two thread-MPI
ranks and four OpenMP threads per rank. You should only use
these options when seeking optimal performance, and
must take care that the ranks you create can have
all of their OpenMP threads run on the same socket.
The number of ranks should be a multiple of the number of
sockets, and the number of cores per node should be
a multiple of the number of threads per rank.

::

    gmx mdrun -ntmpi 4 -nb gpu -pme cpu

Starts :ref:`mdrun <gmx mdrun>` using four thread-MPI ranks. The CPU
cores available will be split evenly between the ranks using OpenMP
threads. The long-range component of the forces are calculated on
CPUs. This may be optimal on hardware where the CPUs are relatively
powerful compared to the GPUs. The bonded part of force calculation
will automatically be assigned to the GPU, since the long-range
component of the forces are calculated on CPU(s).

::

    gmx mdrun -ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu

Starts :ref:`mdrun <gmx mdrun>` using a single thread-MPI rank that
will use all available CPU cores. All interaction types that can run
on a GPU will do so. This may be optimal on hardware where the CPUs
are extremely weak compared to the GPUs.

::

    gmx mdrun -ntmpi 4 -nb gpu -pme cpu -gputasks 0011

Starts :ref:`mdrun <gmx mdrun>` using four thread-MPI ranks, and maps them
to GPUs with IDs 0 and 1. The CPU cores available will be split evenly between
the ranks using OpenMP threads, with the first two ranks offloading short-range
nonbonded force calculations to GPU 0, and the last two ranks offloading to GPU 1.
The long-range component of the forces are calculated on CPUs. This may be optimal
on hardware where the CPUs are relatively powerful compared to the GPUs.

::

    gmx mdrun -ntmpi 4 -nb gpu -pme gpu -npme 1 -gputasks 0001

Starts :ref:`mdrun <gmx mdrun>` using four thread-MPI ranks, one of which is
dedicated to the long-range PME calculation. The first 3 threads offload their
short-range non-bonded calculations to the GPU with ID 0, the 4th (PME) thread
offloads its calculations to the GPU with ID 1.

::

    gmx mdrun -ntmpi 4 -nb gpu -pme gpu -npme 1 -gputasks 0011

Similar to the above example, with 3 ranks assigned to calculating short-range
non-bonded forces, and one rank assigned to calculate the long-range forces.
In this case, 2 of the 3 short-range ranks offload their nonbonded force
calculations to GPU 0. The GPU with ID 1 calculates the short-ranged forces of
the 3rd short-range rank, as well as the long-range forces of the PME-dedicated
rank. Whether this or the above example is optimal will depend on the capabilities
of the individual GPUs and the system composition.

::

    gmx mdrun -gpu_id 12

Starts :ref:`mdrun <gmx mdrun>` using GPUs with IDs 1 and 2 (e.g. because
GPU 0 is dedicated to running a display). This requires
two thread-MPI ranks, and will split the available
CPU cores between them using OpenMP threads.

::

    gmx mdrun -nt 6 -pin on -pinoffset 0 -pinstride 1
    gmx mdrun -nt 6 -pin on -pinoffset 6 -pinstride 1

Starts two :ref:`mdrun <gmx mdrun>` processes, each with six total threads
arranged so that the processes affect each other as little as possible by
being assigned to disjoint sets of physical cores.
Threads will have their affinities set to particular
logical cores, beginning from the first and 7th logical cores, respectively. The
above would work well on an Intel CPU with six physical cores and
hyper-threading enabled. Use this kind of setup only
if restricting :ref:`mdrun <gmx mdrun>` to a subset of cores to share a
node with other processes.
A word of caution: The mapping of logical CPUs/cores to physical
cores may differ between operating systems. On Linux,
``cat /proc/cpuinfo`` can be examined to determine this mapping.

::

    mpirun -np 2 gmx_mpi mdrun

When using an :ref:`gmx mdrun` compiled with external MPI,
this will start two ranks and as many OpenMP threads
as the hardware and MPI setup will permit. If the
MPI setup is restricted to one node, then the resulting
:ref:`gmx mdrun` will be local to that node.

.. _gmx-mdrun-multiple-nodes:

Running :ref:`mdrun <gmx mdrun>` on more than one node
------------------------------------------------------

This requires configuring |Gromacs| to build with an external MPI
library. By default, this :ref:`mdrun <gmx mdrun>` executable is run with
``gmx_mpi mdrun``. All of the considerations for running single-node
:ref:`mdrun <gmx mdrun>` still apply, except that ``-ntmpi`` and ``-nt`` cause a fatal
error, and instead the number of ranks is controlled by the
MPI environment.
Settings such as ``-npme`` are much more important when
using multiple nodes. Configuring the MPI environment to
produce one rank per core is generally good until one
approaches the strong-scaling limit. At that point, using
OpenMP to spread the work of an MPI rank over more than one
core is needed to continue to improve absolute performance.
The location of the scaling limit depends on the processor,
presence of GPUs, network, and simulation algorithm, but
it is worth measuring at around ~200 particles/core if you
need maximum throughput.

There are further command-line parameters that are relevant in these
cases.

``-tunepme``
    Defaults to "on." If "on," a simulation will
    optimize various aspects of the PME and DD algorithms, shifting
    load between ranks and/or GPUs to maximize throughput. Some
    :ref:`mdrun <gmx mdrun>` features are not compatible with this, and these ignore
    this option.

``-dlb``
    Can be set to "auto," "no," or "yes."
    Defaults to "auto." Doing Dynamic Load Balancing between MPI ranks
    is needed to maximize performance. This is particularly important
    for molecular systems with heterogeneous particle or interaction
    density. When a certain threshold for performance loss is
    exceeded, DLB activates and shifts particles between ranks to improve
    performance. If available, using ``-bonded gpu`` is expected
    to improve the ability of DLB to maximize performance.
    DLB is not compatible with GPU-resident parallelization (with ``-update gpu``)
    and therefore it remains switched off in such simulations.

During the simulation :ref:`gmx mdrun` must communicate between all
PP ranks to compute quantities such as kinetic energy for log file
reporting, or perhaps temperature coupling. By default, this happens
whenever necessary to honor several :ref:`mdp options <mdp-general>`,
so that the period between communication phases is the least common
denominator of :mdp:`nstcalcenergy`,
:mdp:`nsttcouple`, and :mdp:`nstpcouple`.

Note that ``-tunepme`` has more effect when there is more than one
:term:`node`, because the cost of communication for the PP and PME
ranks differs. It still shifts load between PP and PME ranks, but does
not change the number of separate PME ranks in use.

Note also that ``-dlb`` and ``-tunepme`` can interfere with each other, so
if you experience performance variation that could result from this,
you may wish to tune PME separately, and run the result with ``mdrun
-notunepme -dlb yes``.

The :ref:`gmx tune_pme` utility is available to search a wider
range of parameter space, including making safe
modifications to the :ref:`tpr` file, and varying ``-npme``.
It is only aware of the number of ranks created by
the MPI environment, and does not explicitly manage
any aspect of OpenMP during the optimization.

Examples for :ref:`mdrun <gmx mdrun>` on more than one node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The examples and explanations for for single-node :ref:`mdrun <gmx mdrun>` are
still relevant, but ``-ntmpi`` is no longer the way
to choose the number of MPI ranks.

::

    mpirun -np 16 gmx_mpi mdrun

Starts :ref:`gmx mdrun` with 16 ranks, which are mapped to
the hardware by the MPI library, e.g. as specified
in an MPI hostfile. The available cores will be
automatically split among ranks using OpenMP threads,
depending on the hardware and any environment settings
such as ``OMP_NUM_THREADS``.

::

    mpirun -np 16 gmx_mpi mdrun -npme 5

Starts :ref:`gmx mdrun` with 16 ranks, as above, and
require that 5 of them are dedicated to the PME
component.

::

    mpirun -np 11 gmx_mpi mdrun -ntomp 2 -npme 6 -ntomp_pme 1

Starts :ref:`gmx mdrun` with 11 ranks, as above, and
require that six of them are dedicated to the PME
component with one OpenMP thread each. The remaining
five do the PP component, with two OpenMP threads
each.

::

    mpirun -np 4 gmx_mpi mdrun -ntomp 6 -nb gpu -gputasks 00

Starts :ref:`gmx mdrun` on a machine with two nodes, using
four total ranks, each rank with six OpenMP threads,
and both ranks on a node sharing GPU with ID 0.

::

    mpirun -np 8 gmx_mpi mdrun -ntomp 3 -gputasks 0000

Using a same/similar hardware as above,
starts :ref:`gmx mdrun` on a machine with two nodes, using
eight total ranks, each rank with three OpenMP threads,
and all four ranks on a node sharing GPU with ID 0.
This may or may not be faster than the previous setup
on the same hardware.

::

    mpirun -np 20 gmx_mpi mdrun -ntomp 4 -gputasks 00

Starts :ref:`gmx mdrun` with 20 ranks, and assigns the CPU cores evenly
across ranks each to one OpenMP thread. This setup is likely to be
suitable when there are ten nodes, each with one GPU, and each node
has two sockets each of four cores.

::

    mpirun -np 10 gmx_mpi mdrun -gpu_id 1

Starts :ref:`gmx mdrun` with 20 ranks, and assigns the CPU cores evenly
across ranks each to one OpenMP thread. This setup is likely to be
suitable when there are ten nodes, each with two GPUs, but another
job on each node is using GPU 0. The job scheduler should set the
affinity of threads of both jobs to their allocated cores, or the
performance of :ref:`mdrun <gmx mdrun>` will suffer greatly.

::

    mpirun -np 20 gmx_mpi mdrun -gpu_id 01

Starts :ref:`gmx mdrun` with 20 ranks. This setup is likely
to be suitable when there are ten nodes, each with two
GPUs, but there is no need to specify ``-gpu_id`` for the
normal case where all the GPUs on the node are available
for use.

Avoiding communication for constraints
--------------------------------------

Because of the very short time it takes to perform an MD step,
in particular close to the scaling limit, any communication will
have a negative effect on performance due to latency overhead
and synchronization. Most of the communication can not be avoided,
but sometimes one can completely avoid communication of coordinates
for constraints. The points listed below will improve performance
in general and can have a particularly strong effect at the scaling
limit which is around ~100 atoms/core or ~10000 atoms/GPU. Simulations
that need to be done as fast as possible, or strong-scaling benchmarks
should be constructed with these points in mind.

When possible, one should avoid the use of ``constraints = all-bonds``
with P-LINCS. This not only requires a lot of communication, it also
sets an artificial minimum on the size of domains. If you are using
an atomistic force field and integrating with a time step of 2 fs,
you can usually change to constraints ``constraints = h-bonds``
without changing other settings. These are
actually the settings most force fields were parameterized with,
so this is also scientifically better.

To completely avoid communication for constraints and/or to have
the update run on a GPU, the system needs to support so-called
"update groups" (or no constraints at all). Update groups are
supported when all atoms involved in coupled constraints are
coupled directly to one central atom and consecutively ordered,
not interdispersed with non-constrained atoms. An example is a
compactly described methyl group. For atomistic
force fields with ``constraints = h-bonds`` this means in practice
that in the topology hydrogens come adjacent to their connected heavy atom.
In addition, when virtual sites are present,
the constructing atoms should all be constrained together and
the virtual site and constructing atoms should be consecutive,
but the order does not matter.
The TIP4P water model is an example of this.
Whether or not update groups are used is noted in the log file.
When they cannot be used, the reason for disabling them is also noted.

Finding out how to run :ref:`mdrun <gmx mdrun>` better
------------------------------------------------------

The Wallcycle module is used for runtime performance measurement of :ref:`gmx mdrun`.
At the end of the log file of each run, the "Real cycle and time accounting" section
provides a table with runtime statistics for different parts of the :ref:`gmx mdrun` code
in rows of the table.
The table contains columns indicating the number of ranks and threads that
executed the respective part of the run, wall-time and cycle
count aggregates (across all threads and ranks) averaged over the entire run.
The last column also shows what percentage of the total runtime each row represents.
Note that the :ref:`gmx mdrun` timer resetting functionalities (``-resethway`` and ``-resetstep``)
reset the performance counters and therefore are useful to avoid startup overhead and
performance instability (e.g. due to load balancing) at the beginning of the run.

The performance counters are:

* Particle-particle during Particle mesh Ewald
* Domain decomposition
* Domain decomposition communication load
* Domain decomposition communication bounds
* Virtual site constraints
* Send X to Particle mesh Ewald
* Neighbor search
* Launch GPU operations
* Communication of coordinates
* Force
* Waiting + Communication of force
* Particle mesh Ewald
* PME redist. X/F
* PME spread
* PME gather
* PME 3D-FFT
* PME 3D-FFT Communication
* PME solve Lennard-Jones
* PME solve LJ
* PME solve Elec
* PME wait for particle-particle
* Wait + Receive PME force
* Wait GPU nonlocal
* Wait GPU local
* Wait PME GPU spread
* Wait PME GPU gather
* Reduce PME GPU Force
* Non-bonded position/force buffer operations
* Virtual site spread
* COM pull force
* AWH (accelerated weight histogram method)
* Write trajectory
* Update
* Constraints
* Communication of energies
* Enforced rotation
* Add rotational forces
* Position swapping
* Interactive MD
* MD Graph

As performance data is collected for every run, they are essential to assessing
and tuning the performance of :ref:`gmx mdrun` performance. Therefore, they benefit
both code developers as well as users of the program.
The counters are an average of the time/cycles different parts of the simulation take,
hence can not directly reveal fluctuations during a single run (although comparisons across
multiple runs are still very useful).

Counters will appear in an MD log file only if the related parts of the code were
executed during the :ref:`gmx mdrun` run. There is also a special counter called "Rest" which
indicates the amount of time not accounted for by any of the counters above. Therefore,
a significant amount "Rest" time (more than a few percent) will often be an indication of
parallelization inefficiency (e.g. serial code) and it is recommended to be reported to the
developers.

An additional set of subcounters can offer more fine-grained inspection of performance. They are:

* Domain decomposition redistribution
* DD neighbor search grid + sort
* DD setup communication
* DD make topology
* DD make constraints
* DD topology other
* Neighbor search grid local
* NS grid non-local
* NS search local
* NS search non-local
* Bonded force
* Bonded-FEP force
* Restraints force
* Listed buffer operations
* Nonbonded pruning
* Nonbonded force
* Launch non-bonded GPU tasks
* Launch PME GPU tasks
* Ewald force correction
* Non-bonded position buffer operations
* Non-bonded force buffer operations

Subcounters are geared toward developers and have to be enabled during compilation. See
:doc:`/dev-manual/build-system` for more information.

..  todo::

    In future patch:
    - red flags in log files, how to interpret wallcycle output
    - hints to devs how to extend wallcycles

.. _gmx-mdrun-on-gpu:

Running :ref:`mdrun <gmx mdrun>` with GPUs
------------------------------------------

.. _gmx-gpu-tasks:

Types of GPU tasks
^^^^^^^^^^^^^^^^^^

To better understand the later sections on different GPU use cases for
calculation of :ref:`short range<gmx-gpu-pp>`, :ref:`PME<gmx-gpu-pme>`,
:ref:`bonded interactions<gmx-gpu-bonded>` and
:ref:`update and constraints <gmx-gpu-update>`
we first introduce the concept of different GPU tasks. When thinking about
running a simulation, several different kinds of interactions between the atoms
have to be calculated (for more information please refer to the reference manual).
The calculation can thus be split into several distinct parts that are largely independent
of each other (hence can be calculated in any order, e.g. sequentially or concurrently),
with the information from each of them combined at the end of
time step to obtain the final forces on each atom and to propagate the system
to the next time point. For a better understanding also please see the section
on :ref:`domain decomposition <gmx-domain-decomp>`.

Of all calculations required for an MD step,
|Gromacs| aims to optimize performance bottom-up for each step
from the lowest level (SIMD unit, cores, sockets, accelerators, etc.).
Therefore many of the individual computation units are
highly tuned for the lowest level of hardware parallelism: the SIMD units.
Additionally, with GPU accelerators used as *co-processors*, some of the work
can be *offloaded*, that is calculated simultaneously/concurrently with the CPU
on the accelerator device, with the result being communicated to the CPU.
Right now, |Gromacs| supports GPU accelerator offload of two tasks:
the short-range :ref:`nonbonded interactions in real space <gmx-gpu-pp>`,
and :ref:`PME <gmx-gpu-pme>`.

|Gromacs| supports two major offload modes: force-offload and GPU-resident.
The former involves offloading some of or all interaction calculations with integration
on the CPU (hence requiring per-step data movement). In the GPU-resident mode
by offloading integration and constraints (when used) less data movement is
necessary.

The force-offload mode is the more broadly supported GPU-acceleration mode
with short-range nonbonded offload supported on a wide range of GPU accelerators
(NVIDIA, AMD, and Intel). This is compatible with the grand majority of
the features and parallelization modes and can be used to scale to large machines.
Simultaneously offloading both short-range nonbonded and long-range
PME work to GPU accelerators has some restrictions in terms of feature and parallelization
compatibility (please see the :ref:`section below <gmx-pme-gpu-limitations>`).
Offloading (most types of) bonded interactions is supported in CUDA and SYCL.
The GPU-resident mode is supported with CUDA and SYCL, but it has additional limitations as
described in :ref:`the GPU update section <gmx-gpu-update>`.


.. _gmx-gpu-pp:

GPU computation of short range nonbonded interactions
.....................................................

.. todo:: make this more elaborate and include figures

Using the GPU for the short-ranged nonbonded interactions provides
the majority of the available speed-up compared to run using only the CPU.
Here, the GPU acts as an accelerator that can effectively parallelize
this problem and thus reduce the calculation time.

.. _gmx-gpu-pme:

GPU accelerated calculation of PME
..................................

.. todo:: again, extend this and add some actual useful information concerning performance etc...

|Gromacs| allows offloading of the PME calculation
to the GPU, to further reduce the load on the CPU and improve usage overlap between
CPU and GPU. Here, the solving of PME will be performed in addition to the calculation
of the short range interactions on the same GPU as the short range interactions.

.. _gmx-pme-gpu-limitations:

Known limitations
.................

**Please note again the limitations outlined below!**

- Only a PME order of 4 is supported on GPUs.

- Multiple ranks (hence multiple GPUs) computing PME have limited support:
  experimental PME decomposition in hybrid mode (``-pmefft cpu``) with
  CUDA from the 2022 release and full GPU PME decomposition since the
  2023 release with CUDA or SYCL (when |Gromacs| is built with
  :ref:`cuFFTMp <cufftmp installation>` or
  :ref:`HeFFTe <heffte installation>`).

- Only dynamical integrators are supported (ie. leap-frog, Velocity Verlet,
  stochastic dynamics)

- LJ PME is not supported on GPUs.

- When |Gromacs| is built with SYCL using oneAPI for AMD/NVIDIA GPUs, only
  hybrid mode (``-pmefft cpu``) is supported. Fully-offloaded PME is supported
  when using oneAPI for Intel GPUs and hipSYCL for AMD/NVIDIA GPUs.

.. _gmx-gpu-bonded:

GPU accelerated calculation of bonded interactions (CUDA and SYCL)
.......................................................................

.. todo:: again, extend this and add some actual useful information concerning performance etc...

|Gromacs| allows the offloading of the bonded part of the PP
workload to a compatible GPU. This is treated as part of the PP
work, and requires that the short-ranged non-bonded task also runs on
a GPU. Typically, there is a performance advantage to offloading
bonded interactions in particular when the amount of CPU resources per GPU
is relatively little (either because the CPU is weak or there are few CPU
cores assigned to a GPU in a run) or when there are other computations on the CPU.
A typical case for the latter is free-energy calculations.

.. _gmx-gpu-update:

GPU accelerated calculation of constraints and coordinate update (CUDA and SYCL only)
.....................................................................................

.. TODO again, extend this with information on when is GPU update supported

|Gromacs| makes it possible to also perform the coordinate update and (if requested)
constraint calculation on a GPU.
This parallelization mode is referred to as "GPU-resident" as all force and coordinate
data can remain resident on the GPU for a number of steps (typically between temperature/pressure coupling or
neighbor searching steps).
The GPU-resident mode allows executing all (supported) computation of a simulation step on the GPU. 
This has the benefit that there is less coupling between CPU host and GPU and
on typical MD steps data does not need to be transferred between CPU and GPU
in contrast to the force-offload scheme requires coordinates and forces to be transferred
every step between the CPU and GPU.
The GPU-resident scheme however is still able to carry out part of the computation
on the CPU concurrently with GPU calculation.
This helps supporting the broad range of |Gromacs| features not all of which are 
ported to GPUs. At the same time, it also allows improving performance by making 
use of the otherwise mostly idle CPU. It can often be advantageous to move the bonded 
or PME calculation back to the CPU, but the details of this will depending on the
relative performance if the CPU cores paired in a simulation with a GPU.

GPU-resident mode is enabled by default (when supported) with an automatic
fallback to CPU update when the build configuration or simulation settings
are incompatible with it. 
It is possible to change the default behaviour by setting the
``GMX_FORCE_UPDATE_DEFAULT_CPU`` environment variable. In this
case simulations following the default behavior (ie. ``-update auto``)
will run the update on the CPU.

Using this parallelization mode is typically advantageous in cases where a fast GPU is
used with a slower CPU, in particular if there is only single simulation assigned to a GPU.
However, in typical throughput cases where multiple runs are assigned to each GPU,
offloading everything, especially without moving back some of the work to the CPU
can perform worse than the parallelization mode where only force computation is offloaded.

Assigning tasks to GPUs
.......................

Depending on which tasks should be performed on which hardware, different kinds of
calculations can be combined on the same or different GPUs, according to the information
provided for running :ref:`mdrun <gmx mdrun>`.

It is possible to assign the calculation of the different computational tasks to the same GPU, meaning
that they will share the computational resources on the same device, or to different processing units
that will each perform one task each.

One overview over the possible task assignments is given below:

|Gromacs| version 2018:

  Two different types of assignable GPU accelerated tasks are available, (short-range) nonbonded and PME.
  Each PP rank has a nonbnonded task that can be offloaded to a GPU.
  If there is only one rank with a PME task (including if that rank is a
  PME-only rank), then that task can be offloaded to a GPU. Such a PME
  task can run wholly on the GPU, or have its latter stages run only on the CPU.

  Limitations are that PME on GPU does not support PME domain decomposition,
  so that only one PME task can be offloaded to a single GPU
  assigned to a separate PME rank, while the nonbonded can be decomposed and offloaded to multiple GPUs.

|Gromacs| version 2019:

  No new assignable GPU tasks are available, but any bonded interactions
  may run on the same GPU as the short-ranged interactions for a PP task.
  This can be influenced with the ``-bonded`` flag.

|Gromacs| version 2020:

  Update and constraints can run on the same GPU as the short-ranged nonbonded and bonded interactions for a PP task.
  This can be influenced with the ``-update`` flag.

|Gromacs| version 2021/2022:

  Communication and auxiliary tasks can also be offloaded in CUDA builds.
  In domain-decomposition halo exchange and PP-PME communication,
  instead of staging transfers between GPUs though the CPU,
  direct GPU--GPU communication is possible.
  As an auxiliary tasks for halo exchange  data packing and unpacking is performed 
  which is also offloaded to the GPU.
  In the 2021 release this is supported with thread-MPI and from the 2022 release
  it is also supported using GPU-aware MPI.
  Direct GPU communication is not enabled by default and can be triggered using
  the ``GMX_ENABLE_DIRECT_GPU_COMM`` environment variable (will only have an effect
  on supported systems).

|Gromacs| version 2023:

  Update now runs by default on the GPU with supported simulation settings; note that this is only available with CUDA and SYCL not with OpenCL.
  
  PME decomposition support adds additional parallelization-related auxiliary GPU tasks including grid packing and reduction operations
  as well as distributed GPU FFT computation.

  Experimental support for CUDA-graphs scheduling has been added, which supports most GPU-resident runs that don't require CPU force computation.


Performance considerations for GPU tasks
........................................

#) The performance balance depends on the speed and number of CPU cores you
   have vs the speed and number of GPUs you have.

#) The GPU-resident parallelization mode (with update/constraints offloaded) is less
   sensitive to the appropriate CPU-GPU balance than the force-offload mode. 

#) With slow/old GPUs and/or fast/modern CPUs with many
   cores, it might make more sense to let the CPU do PME calculation,
   with the GPUs focused on the nonbonded calculation.

#) With fast/modern GPUs and/or slow/old CPUs with few cores,
   it generally helps to have the GPU do PME.

#) Offloading bonded work to a GPU will often not improve simulation performance
   as efficient CPU-based kernels can complete the bonded computation
   before the GPU is done with other offloaded work. Therefore,
   `gmx mdrun` will default to no bonded offload when PME is offloaded.
   Typical cases where performance can improve with bonded offload are:
   with significant bonded work (e.g. pure lipid or mostly polymer systems with little solvent),
   with very few and/or slow CPU cores per GPU, or when the CPU does
   other computation (e.g. PME, free energy).

#) On most modern hardware GPU-resident mode (default) is faster than force-offload mode,
   although it may leave the CPU idle. Moving back the bonded work to the CPU (``-bonded cpu``) is a
   better way to make use of a fast CPU than leaving integration and constraints on the CPU.
   The only exception may be multi-simulations with a significant number of simulations assigned to each GPU.

#) Direct GPU communication will in most cases outperform staged communication (both with thread-MPI and MPI).
   Ideally it should be combined with GPU-resident mode to maximize the benefit.

#) The only way to know for sure which alternative is best for
   your machine is to test and check performance.

.. todo:: we need to be more concrete here, i.e. what machine/software aspects to take into consideration, when will default run mode be using PME-GPU and when will it not, when/how should the user reason about testing different settings than the default.

.. todo:: someone who knows about the mixed mode should comment further.

Reducing overheads in GPU accelerated runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order for CPU cores and GPU(s) to execute concurrently, tasks are
launched and executed asynchronously on the GPU(s) while the CPU cores
execute non-offloaded force computation (like computing bonded forces or free energy computation).
Asynchronous task launches are handled by the GPU device driver and
require CPU involvement. Therefore, scheduling
GPU tasks requires CPU resources that can compete with other CPU tasks
and cause interference that could lead to slowdown.

Delays in CPU execution are caused by the latency of launching GPU tasks,
an overhead that can become significant as simulation ns/day increases
(i.e. with shorter wall-time per step).
The cost of launching GPU work is measured by :ref:`gmx mdrun` and reported in the performance
summary section of the log file ("Launch PP GPU ops."/"Launch PME GPU ops." rows).
A few percent of runtime spent in launching work is normal,
but in fast-iterating and multi-GPU parallel runs, costs of 10% or larger can be observed.
Whether this has a significant performance impact depends on how much work
within the main MD step is assigned to the CPU. With most or all force computation offloaded,
and when the CPU is not involved in communication (e.g. with thread-MPI and direct GPU communication enabled) 
it may be that large launch costs do not lead to large performance losses.
However, when the CPU is assigned computation (e.g. in free energy or pull/AWH simulations)
or MPI communication is launched from the CPU (even with GPU-aware MPI), the
GPU launch cost will compete with other CPU work and therefore represent overheads.
In general, a user can do little to avoid such overheads, but there
are a few cases where tweaks can give performance benefits.
In OpenCL runs, timing of GPU tasks is by default enabled and,
while in most cases its impact is small, in fast runs performance can be affected.
In these cases, when more than a few percent of "Launch GPU ops" time is observed,
it is recommended to turn off timing by setting the ``GMX_DISABLE_GPU_TIMING``
environment variable.
In parallel runs with many ranks sharing a GPU,
launch overheads can also be reduced by starting fewer thread-MPI
or MPI ranks per GPU; e.g. most often one rank per thread or core is not optimal.
The CUDA graphs functionality (added in |Gromacs| 2023) targets reducing such
overheads and improving GPU work scheduling efficiency and therefore
it can provide significant improvements especially for small simulation systems
running on fast GPUs. Since it is a new feature, in the 2023 release CUDA-graph support
needs to be triggered using the ``GMX_CUDA_GRAPH`` environment variable.

The second type of overhead, interference of the GPU runtime or driver with CPU computation,
is caused by the scheduling and coordination of GPU tasks.
A separate GPU runtime/driver thread requires CPU resources
which may compete with the concurrently running non-offloaded tasks (if present),
potentially degrading the performance of this CPU work.
To minimize the overhead it can be useful to
leave at least one CPU hardware thread unused when launching :ref:`gmx mdrun`,
especially on CPUs with high core counts and/or simultaneous multithreading enabled.
E.g. on a machine with a 16-core CPU and 32 threads,
try ``gmx mdrun -ntomp 31 -pin on``.
This will leave some CPU resources for the GPU task scheduling
potentially reducing interference with CPU computation.
Note that assigning fewer resources to :ref:`gmx mdrun` CPU computation
involves a tradeoff which, with many CPU cores per GPU, may not be significant,
but in some cases (e.g. with multi-rank MPI runs) it may lead to complex
resource assignment and may outweigh the benefits of reduced GPU scheduling overheads,
so we recommend to test the alternatives before adopting such techniques.

.. todo:: In future patch: any tips not covered above

Running the OpenCL version of mdrun
-----------------------------------

Currently supported hardware architectures are:

- GCN-based and CDNA-based AMD GPUs;
- NVIDIA GPUs prior to Volta;
- Intel iGPUs.

Make sure that you have the latest drivers installed. For AMD GPUs,
the compute-oriented `ROCm <https://rocm.docs.amd.com/en/latest/>`_ stack is recommended;
alternatively, the AMDGPU-PRO stack is also compatible; using the outdated
and unsupported ``fglrx`` proprietary driver and runtime is not recommended (but
for certain older hardware that may be the only way to obtain support).
In addition Mesa version 17.0 or newer with LLVM 4.0 or newer is also supported.
For NVIDIA GPUs, using the proprietary driver is
required as the open source nouveau driver (available in Mesa) does not
provide the OpenCL support.
For Intel integrated GPUs, the `Neo driver <https://github.com/intel/compute-runtime/releases>`_ is
recommended.

The minimum OpenCL version required is |REQUIRED_OPENCL_MIN_VERSION|. See
also the :ref:`known limitations <opencl-known-limitations>`.

Devices from the AMD GCN architectures (all series) are compatible
and regularly tested; NVIDIA Kepler and later (compute capability 3.0)
are known to work, but before doing production runs always make sure that the |Gromacs| tests
pass successfully on the hardware.

The OpenCL GPU kernels are compiled at run time. Hence,
building the OpenCL program can take a few seconds, introducing a slight
delay in the :ref:`gmx mdrun` startup. This is not normally a
problem for long production MD, but you might prefer to do some kinds
of work, e.g. that runs very few steps, on just the CPU (e.g. see ``-nb`` above).

The same ``-gpu_id`` option (or ``GMX_GPU_ID`` environment variable)
used to select CUDA devices, or to define a mapping of GPUs to PP
ranks, is used for OpenCL devices.

Some other :ref:`OpenCL management <opencl-management>` environment
variables may be of interest to developers.

.. _opencl-known-limitations:

Known limitations of the OpenCL support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Limitations in the current OpenCL support of interest to |Gromacs| users:

- Intel integrated GPUs are supported. Intel CPUs and Xeon Phi are not supported.
  Set ``-DGMX_GPU_NB_CLUSTER_SIZE=4`` when compiling |Gromacs| to run on consumer
  Intel GPUs (as opposed to Ponte Vecchio / Data Center Max GPUs).
- Due to blocking behavior of some asynchronous task enqueuing functions
  in the NVIDIA OpenCL runtime, with the affected driver versions there is
  almost no performance gain when using NVIDIA GPUs.
  The issue affects NVIDIA driver versions up to 349 series, but it
  known to be fixed 352 and later driver releases.
- On NVIDIA GPUs the OpenCL kernels achieve much lower performance
  than the equivalent CUDA kernels due to limitations of the NVIDIA OpenCL
  compiler.
- On the NVIDIA Volta and Turing architectures the OpenCL code is known to produce
  incorrect results with driver version up to 440.x (most likely due to compiler issues).
  Runs typically fail on these architectures.

Running SYCL version of mdrun
-----------------------------

Make sure that you have the latest drivers installed and check the :ref:`installation guide <SYCL GPU acceleration>`
for the list of compatible hardware and software and the recommended compile-time options.

Please keep in mind the following environment variables that might be useful:

- When using oneAPI runtime:

  - ``SYCL_CACHE_PERSISTENT=1``: enables caching of GPU kernels, reducing :ref:`gmx mdrun` startup time.

In addition to ``-gpu_id`` option, backend-specific environment variables, like ``SYCL_DEVICE_FILTER``
or ``ROCR_VISIBLE_DEVICES``, could be used to select GPUs.

Performance checklist
---------------------

There are many different aspects that affect the performance of simulations in
|Gromacs|. Most simulations require a lot of computational resources, therefore
it can be worthwhile to optimize the use of those resources. Several issues
mentioned in the list below could lead to a performance difference of a factor
of 2. So it can be useful go through the checklist.

|Gromacs| configuration
^^^^^^^^^^^^^^^^^^^^^^^

* Don't use double precision unless you're absolute sure you need it.
* Compile the FFTW library (yourself) with the correct flags on x86 (in most
  cases, the correct flags are automatically configured).
* On x86, use gcc as the compiler (not icc, pgi or the Cray compiler).
* On POWER, use gcc instead of IBM's xlc.
* Use a new compiler version, especially for gcc (e.g. from version 5 to 6
  the performance of the compiled code improved a lot).
* MPI library: OpenMPI usually has good performance and causes little trouble.
* Make sure your compiler supports OpenMP (some versions of Clang don't).
* If you have GPUs that support either CUDA, OpenCL, or SYCL, use them.

  * Configure with ``-DGMX_GPU=CUDA``, ``-DGMX_GPU=OpenCL``, or ``-DGMX_GPU=SYCL``.
  * For CUDA, use the newest CUDA available for your GPU to take advantage of the
    latest performance enhancements.
  * Use a recent GPU driver.
  * Make sure you use an :ref:`gmx mdrun` with ``GMX_SIMD`` appropriate for the CPU
    architecture; the log file will contain a warning note if suboptimal setting is used.
    However, prefer ``AVX2`` over ``AVX512`` in GPU or highly parallel MPI runs (for more
    information see the :ref:`intra-core parallelization information <intra-core-parallelization>`).
  * If compiling on a cluster head node, make sure that ``GMX_SIMD``
    is appropriate for the compute nodes.

Run setup
^^^^^^^^^

* For an approximately spherical solute, use a rhombic dodecahedron unit cell.
* When using a time-step of <=2.5 fs, use :mdp-value:`constraints=h-bonds`
  (and not :mdp-value:`constraints=all-bonds`), since:

  * this is faster, especially with GPUs;
  * it is necessary to be able to use GPU-resident mode;
  * and most force fields have been parametrized with only bonds involving hydrogens constrained.

* You can often increase the time-step to 4 fs by repartitioning hydrogen
  masses using the ``mass-repartition-factor`` mdp option. This does not
  affect equilibrium distributions, but makes dynamics slightly slower.
* You can increase the time-step to 4 or 5 fs when using virtual interaction
  sites (``gmx pdb2gmx -vsite h``).
* For massively parallel runs with PME, you might need to try different numbers
  of PME ranks (``gmx mdrun -npme ???``) to achieve best performance;
  :ref:`gmx tune_pme` can help automate this search.
* For massively parallel runs (also ``gmx mdrun -multidir``), or with a slow
  network, global communication can become a bottleneck and you can reduce it
  by choosing larger periods for algorithms such as temperature and
  pressure coupling).

Checking and improving performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Look at the end of the ``md.log`` file to see the performance and the cycle
  counters and wall-clock time for different parts of the MD calculation. The
  PP/PME load ratio is also printed, with a warning when a lot of performance is
  lost due to imbalance.
* Adjust the number of PME ranks and/or the cut-off and PME grid-spacing when
  there is a large PP/PME imbalance. Note that even with a small reported
  imbalance, the automated PME-tuning might have reduced the initial imbalance.
  You could still gain performance by changing the mdp parameters or increasing
  the number of PME ranks.
* (Especially) In GPU-resident runs (``-update gpu``):

  * Frequent virial or energy computation can have a large overhead (and this will not show up in the cycle counters).
    To reduce this overhead, increase ``nstcalcenergy``;
  * Frequent temperature or pressure coupling can have significant overhead; 
    to reduce this, make sure to have as infrequent coupling as your algorithms allow (typically >=50-100 steps).

* If the neighbor searching and/or domain decomposition takes a lot of time, increase ``nstlist``. If a Verlet
  buffer tolerance is used, this is done automatically by :ref:`gmx mdrun`
  and the pair-list buffer is increased to keep the energy drift constant.

    * especially with multi-GPU runs, the automatic increasing of ``nstlist`` at ``mdrun``
      startup can be conservative and larger value is often be optimal
      (e.g. ``nstlist=200-300`` with PME and default Verlet buffer tolerance).

    * odd values of nstlist should be avoided when using CUDA Graphs
      to minimize the overhead associated with graph instantiation.

* If ``Comm. energies`` takes a lot of time (a note will be printed in the log
  file), increase ``nstcalcenergy``.
* If all communication takes a lot of time, you might be running on too many
  cores, or you could try running combined MPI/OpenMP parallelization with 2
  or 4 OpenMP threads per MPI process.
* In multi-GPU runs avoid using as many ranks as cores (or hardware threads) since
  this introduces a major inefficiency due to overheads associated to GPUs sharing by several MPI ranks.
  Use at most a few ranks per GPU, 1-3 ranks is generally optimal;
  with GPU-resident mode and direct GPU communication typically 1 rank/GPU is best.

