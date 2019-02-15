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
        the hardware bios or through a setting in the Linux operating
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

    MPI
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
        on AMD devices (both GPUs and APUs) and Intel integrated GPUs; NVIDIA
        hardware is also supported.

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

Core level parallelization via SIMD: SSE, AVX, etc.
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
the ``-DGMX_SIMD`` CMake option can be used. For best performance always pick the highest
(latest) SIMD instruction set supported by the target architecture (and |Gromacs|). To use a single
installation on multiple different machines, it is convenient to compile the analysis tools with
the lowest common SIMD instruction set (as these rely little on SIMD acceleration), but for best
performance :ref:`mdrun <gmx mdrun>` should be compiled separately for each machine.

.. TODO add a note on AVX throttle and its impact on MPI-parallel and GPU accelerated runs

Process(-or) level parallelization via OpenMP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|Gromacs| :ref:`mdrun <gmx mdrun>` supports OpenMP multithreading for all parts
of the code for the :doc:`Verlet cut-off scheme <cutoff-schemes>`, as well as for the PME
code in case of the group scheme. OpenMP is enabled by default and
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

Thread-MPI is included in the |Gromacs| source and it is the default parallelization since
version 4.5, practically rendering the serial :ref:`mdrun <gmx mdrun>` deprecated.
Compilation with thread-MPI is controlled by the ``GMX_THREAD_MPI`` CMake variable.

Thread-MPI is compatible with most :ref:`mdrun <gmx mdrun>` features and parallelization schemes,
including OpenMP, GPUs; it is not compatible with MPI and multi-simulation runs.

By default, the thread-MPI mdrun will use all available cores in the machine by starting
an appropriate number of ranks or OpenMP threads to occupy all of them. The number of
ranks can be controlled using the
``-nt`` and ``-ntmpi`` options. ``-nt`` represents the total number of threads
to be used (which can be a mix of thread-MPI and OpenMP threads with the
:doc:`Verlet scheme <cutoff-schemes>`).

Hybrid/heterogeneous acceleration
.................................

Hybrid acceleration means distributing compute work between available CPUs and GPUs
to improve simulation performance.
Along the :doc:`Verlet cut-off scheme <cutoff-schemes>` new non-bonded algorithms
have been developed with the aim of efficient acceleration both on CPUs and GPUs.

The most compute-intensive parts of simulations, non-bonded force calculation, as well
as possibly the PME and bonded force calculation can be
offloaded to GPUs and carried out simultaneously with remaining CPU work.
Native GPU acceleration is supported with the :doc:`Verlet cut-off scheme <cutoff-schemes>`
(not with the group scheme) with PME, reaction-field, and plain cut-off electrostatics.
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
and there would be no performance gain in the Verlet cut-off scheme.

While the automated CPU-GPU load balancing always attempts to find the optimal cut-off setting,
it might not always be possible to balance CPU and GPU workload. This happens when the CPU threads
finish calculating the bonded forces and PME faster than the GPU the non-bonded force calculation,
even with the shortest possible cut-off. In such cases the CPU will wait for the GPU and this
time will show up as ``Wait GPU local`` in the cycle and timing summary table at the end
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
especially when running separate multi-threaded PME nodes. Depending on the architecture,
input system size, as well as other factors, MPI+OpenMP runs can be as fast and faster
already at small number of processes (e.g. multi-processor Intel Westmere or Sandy Bridge),
but can also be considerably slower (e.g. multi-processor AMD Interlagos machines). However,
there is a more pronounced benefit of multi-level parallelization in highly parallel runs.

Separate PME nodes
^^^^^^^^^^^^^^^^^^

On CPU nodes, particle-particle (PP) and PME calculations are done in the same process one after
another. As PME requires all-to-all global communication, this is most of the time the limiting
factor to scaling on a large number of cores. By designating a subset of nodes for PME
calculations only, performance of parallel runs can be greatly improved.

OpenMP mutithreading in PME nodes is also possible and is supported with both group and
Verlet cut-off schemes. Using multi-threading in PME can can improve performance at high
parallelization. The reason for this is that with N>1 threads the number of processes
communicating, and therefore the number of messages, is reduced by a factor of N.
But note that modern communication networks can process several messages simultaneously,
such that it could be advantageous to have more processes communicating.

Separate PME nodes are not used at low parallelization, the switch at higher parallelization
happens automatically (at > 16 processes). The number of PME nodes is estimated by mdrun.
If the PME load is higher than the PP load, mdrun will automatically balance the load, but
this leads to additional (non-bonded) calculations. This avoids the idling of a large fraction
of the nodes; usually 3/4 of the nodes are PP nodes. But to ensure the best absolute performance
of highly parallel runs, it is advisable to tweak this number which is automated by
the :ref:`tune_pme <gmx tune_pme>` tool.

The number of PME nodes can be set manually on the :ref:`mdrun <gmx mdrun>` command line using the ``-npme``
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
    Setting "gpu" requires that a compatible GPU is available and will be used.
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

    gmx mdrun -ntmpi 1 -nb gpu -pme gpu -bonded gpu

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
:ref:`mdrun_mpi`. All of the considerations for running single-node
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
    Defaults to "on." If "on," a Verlet-scheme simulation will
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

``-gcom``
    During the simulation :ref:`gmx mdrun` must communicate between all ranks to
    compute quantities such as kinetic energy. By default, this
    happens whenever plausible, and is influenced by a lot of
    :ref:`mdp options. <mdp-general>` The period between communication phases
    must be a multiple of :mdp:`nstlist`, and defaults to
    the minimum of :mdp:`nstcalcenergy` and :mdp:`nstlist`.
    ``mdrun -gcom`` sets the number of steps that must elapse between
    such communication phases, which can improve performance when
    running on a lot of ranks. Note that this means that _e.g._
    temperature coupling algorithms will
    effectively remain at constant energy until the next
    communication phase. :ref:`gmx mdrun` will always honor the
    setting of ``mdrun -gcom``, by changing :mdp:`nstcalcenergy`,
    :mdp:`nstenergy`, :mdp:`nstlog`, :mdp:`nsttcouple` and/or
    :mdp:`nstpcouple` if necessary.

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

Starts :ref:`mdrun_mpi` with 16 ranks, which are mapped to
the hardware by the MPI library, e.g. as specified
in an MPI hostfile. The available cores will be
automatically split among ranks using OpenMP threads,
depending on the hardware and any environment settings
such as ``OMP_NUM_THREADS``.

::

    mpirun -np 16 gmx_mpi mdrun -npme 5

Starts :ref:`mdrun_mpi` with 16 ranks, as above, and
require that 5 of them are dedicated to the PME
component.

::

    mpirun -np 11 gmx_mpi mdrun -ntomp 2 -npme 6 -ntomp_pme 1

Starts :ref:`mdrun_mpi` with 11 ranks, as above, and
require that six of them are dedicated to the PME
component with one OpenMP thread each. The remaining
five do the PP component, with two OpenMP threads
each.

::

    mpirun -np 4 gmx_mpi mdrun -ntomp 6 -nb gpu -gputasks 00

Starts :ref:`mdrun_mpi` on a machine with two nodes, using
four total ranks, each rank with six OpenMP threads,
and both ranks on a node sharing GPU with ID 0.

::

    mpirun -np 8 gmx_mpi mdrun -ntomp 3 -gputasks 0000

Using a same/similar hardware as above,
starts :ref:`mdrun_mpi` on a machine with two nodes, using
eight total ranks, each rank with three OpenMP threads,
and all four ranks on a node sharing GPU with ID 0.
This may or may not be faster than the previous setup
on the same hardware.

::

    mpirun -np 20 gmx_mpi mdrun -ntomp 4 -gputasks 00

Starts :ref:`mdrun_mpi` with 20 ranks, and assigns the CPU cores evenly
across ranks each to one OpenMP thread. This setup is likely to be
suitable when there are ten nodes, each with one GPU, and each node
has two sockets each of four cores.

::

    mpirun -np 10 gmx_mpi mdrun -gpu_id 1

Starts :ref:`mdrun_mpi` with 20 ranks, and assigns the CPU cores evenly
across ranks each to one OpenMP thread. This setup is likely to be
suitable when there are ten nodes, each with two GPUs, but another
job on each node is using GPU 0. The job scheduler should set the
affinity of threads of both jobs to their allocated cores, or the
performance of :ref:`mdrun <gmx mdrun>` will suffer greatly.

::

    mpirun -np 20 gmx_mpi mdrun -gpu_id 01

Starts :ref:`mdrun_mpi` with 20 ranks. This setup is likely
to be suitable when there are ten nodes, each with two
GPUs, but there is no need to specify ``-gpu_id`` for the
normal case where all the GPUs on the node are available
for use.

Approaching the scaling limit
-----------------------------

There are several aspects of running a |Gromacs| simulation that are important as the number
of atoms per core approaches the current scaling limit of ~100 atoms/core.

One of these is that the use of ``constraints = all-bonds``  with P-LINCS
sets an artificial minimum on the size of domains. You should reconsider the use
of constraints to all bonds (and bear in mind possible consequences on the safe maximum for dt),
or change lincs_order and lincs_iter suitably.

Finding out how to run :ref:`mdrun <gmx mdrun>` better
------------------------------------------------------

The Wallcycle module is used for runtime performance measurement of :ref:`gmx mdrun`.
At the end of the log file of each run, the "Real cycle and time accounting" section
provides a table with runtime statistics for different parts of the :ref:`gmx mdrun` code
in rows of the table.
The table contains colums indicating the number of ranks and threads that
executed the respective part of the run, wall-time and cycle
count aggregates (across all threads and ranks) averaged over the entire run.
The last column also shows what precentage of the total runtime each row represents.
Note that the :ref:`gmx mdrun` timer resetting functionalities (`-resethway` and `-resetstep`)
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

.. TODO In future patch:
   - red flags in log files, how to interpret wallcycle output
   - hints to devs how to extend wallcycles

.. _gmx-mdrun-on-gpu:

Running :ref:`mdrun <gmx mdrun>` with GPUs
------------------------------------------

.. _gmx-gpu-tasks:

Types of GPU tasks
^^^^^^^^^^^^^^^^^^

To better understand the later sections on different GPU use cases for
calculation of :ref:`short range<gmx-gpu-pp>` and :ref:`PME <gmx-gpu-pme>`,
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
GROMACS aims to optimize performance bottom-up for each step
from the lowest level (SIMD unit, cores, sockets, accelerators, etc.).
Therefore many of the individual computation units are
highly tuned for the lowest level of hardware parallelism: the SIMD units.
Additionally, with GPU accelerators used as *co-processors*, some of the work
can be *offloaded*, that is calculated simultaneously/concurrently with the CPU
on the accelerator device, with the result being communicated to the CPU.
Right now, |Gromacs| supports GPU accelerator offload of two tasks:
the short-range :ref:`nonbonded interactions in real space <gmx-gpu-pp>`,
and :ref:`PME <gmx-gpu-pme>`.

**Please note that the solving of PME on GPU is still only the initial
version supporting this behaviour, and comes with a set of limitations
outlined further below.**

Right now, we generally support short-range nonbonded offload with and
without dynamic pruning on a wide range of GPU accelerators
(both NVIDIA and AMD). This is compatible with the grand majority of
the features and parallelization modes and can be used to scale to large machines.

Simultaneously offloading both short-range nonbonded and long-range
PME work to GPU accelerators is a new feature that that has some
restrictions in terms of feature and parallelization
compatibility (please see the :ref:`section below <gmx-pme-gpu-limitations>`).

.. _gmx-gpu-pp:

GPU computation of short range nonbonded interactions
.....................................................

.. TODO make this more elaborate and include figures

Using the GPU for the short-ranged nonbonded interactions provides
the majority of the available speed-up compared to run using only the CPU.
Here, the GPU acts as an accelerator that can effectively parallelize
this problem and thus reduce the calculation time.

.. _gmx-gpu-pme:

GPU accelerated calculation of PME
..................................

.. TODO again, extend this and add some actual useful information concerning performance etc...

|Gromacs| now allows the offloading of the PME calculation
to the GPU, to further reduce the load on the CPU and improve usage overlap between
CPU and GPU. Here, the solving of PME will be performed in addition to the calculation
of the short range interactions on the same GPU as the short range interactions.

.. _gmx-pme-gpu-limitations:

Known limitations
.................

**Please note again the limitations outlined below!**

- Only compilation with CUDA is supported.

- Only a PME order of 4 is supported on GPUs.

- PME will run on a GPU only when exactly one rank has a
  PME task, ie. decompositions with multiple ranks doing PME are not supported.

- Only single precision is supported.

- Free energy calculations where charges are perturbed are not supported,
  because only single PME grids can be calculated.

- Only dynamical integrators are supported (ie. leap-frog, Velocity Verlet,
  stochastic dynamics)

- LJ PME is not supported on GPUs.

GPU accelerated calculation of bonded interactions (CUDA only)
..............................................................

.. TODO again, extend this and add some actual useful information concerning performance etc...

|Gromacs| now allows the offloading of the bonded part of the PP
workload to a CUDA-compatible GPU. This is treated as part of the PP
work, and requires that the short-ranged non-bonded task also runs on
a GPU. It is an advantage usually only when the CPU is relatively weak
compared with the GPU, perhaps because its workload is too large for
the available cores. This would likely be the case for free-energy
calculations.

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

  Two different types of assignable GPU accelerated tasks are available, NB and PME.
  Each PP rank has a NB task that can be offloaded to a GPU.
  If there is only one rank with a PME task (including if that rank is a
  PME-only rank), then that task can be offloaded to a GPU. Such a PME
  task can run wholly on the GPU, or have its latter stages run only on the CPU.

  Limitations are that PME on GPU does not support PME domain decomposition,
  so that only one PME task can be offloaded to a single GPU
  assigned to a separate PME rank, while NB can be decomposed and offloaded to multiple GPUs.

|Gromacs| version 2019:

  No new assignable GPU tasks are available, but any bonded interactions
  may run on the same GPU as the short-ranged interactions for a PP task.
  This can be influenced with the ``-bonded`` flag.

Performance considerations for GPU tasks
........................................

#) The performance balance depends on the speed and number of CPU cores you
   have vs the speed and number of GPUs you have.

#) With slow/old GPUs and/or fast/modern CPUs with many
   cores, it might make more sense to let the CPU do PME calculation,
   with the GPUs focused on the calculation of the NB.

#) With fast/modern GPUs and/or slow/old CPUs with few cores,
   it generally helps to have the GPU do PME.

#) Offloading bonded work to a GPU will often not improve simulation performance
   as efficient CPU-based kernels can complete the bonded computation
   before the GPU is done with other offloaded work. Therefore,
   `gmx mdrun` will default to no bonded offload when PME is offloaded.
   Typical cases where performance can be improvement with bonded offload are:
   with significant bonded work (e.g. pure lipid or mostly polymer systems with little solvent),
   with very few and/or slow CPU cores per GPU, or when the CPU does
   other computation (e.g. PME, free energy).

#) It *is* possible to use multiple GPUs with PME offload
   by letting e.g.
   3 MPI ranks use one GPU each for short-range interactions,
   while a fourth rank does the PME on its GPU.

#) The only way to know for sure what alternative is best for
   your machine is to test and check performance.

.. TODO: we need to be more concrete here, i.e. what machine/software aspects to take into consideration, when will default run mode be using PME-GPU and when will it not, when/how should the user reason about testing different settings than the default.

.. TODO someone who knows about the mixed mode should comment further.

Reducing overheads in GPU accelerated runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order for CPU cores and GPU(s) to execute concurrently, tasks are
launched and executed asynchronously on the GPU(s) while the CPU cores
execute non-offloaded force computation (like long-range PME electrostatics).
Asynchronous task launches are handled by GPU device driver and
require CPU involvement. Therefore, the work of scheduling
GPU tasks will incur an overhead that can in some cases significantly
delay or interfere with the CPU execution.

Delays in CPU execution are caused by the latency of launching GPU tasks,
an overhead that can become significant as simulation ns/day increases
(i.e. with shorter wall-time per step).
The overhead is measured by :ref:`gmx mdrun` and reported in the performance
summary section of the log file ("Launch GPU ops" row).
A few percent of runtime spent in this category is normal,
but in fast-iterating and multi-GPU parallel runs 10% or larger overheads can be observed.
In general, a user can do little to avoid such overheads, but there
are a few cases where tweaks can give performance benefits.
In single-rank runs timing of GPU tasks is by default enabled and,
while in most cases its impact is small, in fast runs performance can be affected.
The performance impact will be most significant on NVIDIA GPUs with CUDA,
less on AMD and Intel with OpenCL.
In these cases, when more than a few percent of "Launch GPU ops" time is observed,
it is recommended to turn off timing by setting the ``GMX_DISABLE_GPU_TIMING``
environment variable.
In parallel runs with many ranks sharing a GPU,
launch overheads can also be reduced by starting fewer thread-MPI
or MPI ranks per GPU; e.g. most often one rank per thread or core is not optimal.

The second type of overhead, interference of the GPU driver with CPU computation,
is caused by the scheduling and coordination of GPU tasks.
A separate GPU driver thread can require CPU resources
which may clash with the concurrently running non-offloaded tasks,
potentially degrading the performance of PME or bonded force computation.
This effect is most pronounced when using AMD GPUs with OpenCL with
older driver releases (e.g. fglrx 12.15).
To minimize the overhead it is recommended to
leave a CPU hardware thread unused when launching :ref:`gmx mdrun`,
especially on CPUs with high core counts and/or HyperThreading enabled.
E.g. on a machine with a 4-core CPU and eight threads (via HyperThreading) and an AMD GPU,
try ``gmx mdrun -ntomp 7 -pin on``.
This will leave free CPU resources for the GPU task scheduling
reducing interference with CPU computation.
Note that assigning fewer resources to :ref:`gmx mdrun` CPU computation
involves a tradeoff which may outweigh the benefits of reduced GPU driver overhead,
in particular without HyperThreading and with few CPU cores.

.. TODO In future patch: any tips not covered above

Running the OpenCL version of mdrun
-----------------------------------

Currently supported hardware architectures are:
- GCN-based AMD GPUs;
- NVIDIA GPUs (with at least OpenCL 1.2 support);
- Intel iGPUs.
Make sure that you have the latest drivers installed. For AMD GPUs,
the compute-oriented `ROCm <https://rocm.github.io/>`_ stack is recommended;
alternatively, the AMDGPU-PRO stack is also compatible; using the outdated
and unsupported `fglrx` proprietary driver and runtime is not recommended (but
for certain older hardware that may be the only way to obtain support).
In addition Mesa version 17.0 or newer with LLVM 4.0 or newer is also supported.
For NVIDIA GPUs, using the proprietary driver is
required as the open source nouveau driver (available in Mesa) does not
provide the OpenCL support.
For Intel integrated GPUs, the `Neo driver <https://github.com/intel/compute-runtime/releases>`_ is
recommended.
TODO: add more Intel driver recommendations
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
- Due to blocking behavior of some asynchronous task enqueuing functions
  in the NVIDIA OpenCL runtime, with the affected driver versions there is
  almost no performance gain when using NVIDIA GPUs.
  The issue affects NVIDIA driver versions up to 349 series, but it
  known to be fixed 352 and later driver releases.
- On NVIDIA GPUs the OpenCL kernels achieve much lower performance
  than the equivalent CUDA kernels due to limitations of the NVIDIA OpenCL
  compiler.
- PME is currently only supported on AMD devices, because of known
  issues with devices from other vendors

Limitations of interest to |Gromacs| developers:

- The current implementation is not compatible with OpenCL devices that are
  not using warp/wavefronts or for which the warp/wavefront size is not a
  multiple of 32

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
* On x86, use gcc or icc as the compiler (not pgi or the Cray compiler).
* On POWER, use gcc instead of IBM's xlc.
* Use a new compiler version, especially for gcc (e.g. from version 5 to 6
  the performance of the compiled code improved a lot).
* MPI library: OpenMPI usually has good performance and causes little trouble.
* Make sure your compiler supports OpenMP (some versions of Clang don't).
* If you have GPUs that support either CUDA or OpenCL, use them.

  * Configure with ``-DGMX_GPU=ON`` (add ``-DGMX_USE_OPENCL=ON`` for OpenCL).
  * For CUDA, use the newest CUDA availabe for your GPU to take advantage of the
    latest performance enhancements.
  * Use a recent GPU driver.
  * If compiling on a cluster head node, make sure that ``GMX_SIMD``
    is appropriate for the compute nodes.

Run setup
^^^^^^^^^

* For an approximately spherical solute, use a rhombic dodecahedron unit cell.
* When using a time-step of 2 fs, use :mdp-value:`constraints=h-bonds`
  (and not :mdp-value:`constraints=all-bonds`), since this is faster, especially with GPUs,
  and most force fields have been parametrized with only bonds involving
  hydrogens constrained.
* You can increase the time-step to 4 or 5 fs when using virtual interaction
  sites (``gmx pdb2gmx -vsite h``).
* For massively parallel runs with PME, you might need to try different numbers
  of PME ranks (``gmx mdrun -npme ???``) to achieve best performance;
  :ref:`gmx tune_pme` can help automate this search.
* For massively parallel runs (also ``gmx mdrun -multidir``), or with a slow
  network, global communication can become a bottleneck and you can reduce it
  with ``gmx mdrun -gcom`` (note that this does affect the frequency of
  temperature and pressure coupling).

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
* If the neighbor searching takes a lot of time, increase nstlist (with the
  Verlet cut-off scheme, this automatically adjusts the size of the neighbour
  list to do more non-bonded computation to keep energy drift constant).

  * If ``Comm. energies`` takes a lot of time (a note will be printed in the log
    file), increase nstcalcenergy or use ``mdrun -gcom``.
  * If all communication takes a lot of time, you might be running on too many
    cores, or you could try running combined MPI/OpenMP parallelization with 2
    or 4 OpenMP threads per MPI process.
