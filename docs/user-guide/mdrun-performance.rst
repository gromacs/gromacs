Getting good performance from mdrun
===================================
The |Gromacs| build system and the :ref:`gmx mdrun` tool has a lot of built-in
and configurable intelligence to detect your hardware and make pretty
effective use of that hardware. For a lot of casual and serious use of
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
        on AMD devices (both GPUs and APUs); NVIDIA hardware is also supported.

    SIMD
        Modern CPU cores have instructions that can execute large
        numbers of floating-point instructions in a single cycle.


|Gromacs| background information
--------------------------------
The algorithms in :ref:`gmx mdrun` and their implementations are most relevant
when choosing how to make good use of the hardware. For details,
see the Reference Manual. The most important of these are

.. glossary::

    Domain Decomposition
        The domain decomposition (DD) algorithm decomposes the
        (short-ranged) component of the non-bonded interactions into
        domains that share spatial locality, which permits the use of
        efficient algorithms. Each domain handles all of the
        particle-particle (PP) interactions for its members, and is
        mapped to a single MPI rank. Within a PP rank, OpenMP threads
        can share the workload, and some work can be off-loaded to a
        GPU. The PP rank also handles any bonded interactions for the
        members of its domain. A GPU may perform work for more than
        one PP rank, but it is normally most efficient to use a single
        PP rank per GPU and for that rank to have thousands of
        particles. When the work of a PP rank is done on the CPU, mdrun
        will make extensive use of the SIMD capabilities of the
        core. There are various `command-line options
        <controlling-the-domain-decomposition-algorithm` to control
        the behaviour of the DD algorithm.

    Particle-mesh Ewald
        The particle-mesh Ewald (PME) algorithm treats the long-ranged
        components of the non-bonded interactions (Coulomb and/or
        Lennard-Jones).  Either all, or just a subset of ranks may
        participate in the work for computing long-ranged component
        (often inaccurately called simple the "PME"
        component). Because the algorithm uses a 3D FFT that requires
        global communication, its performance gets worse as more ranks
        participate, which can mean it is fastest to use just a subset
        of ranks (e.g.  one-quarter to one-half of the ranks). If
        there are separate PME ranks, then the remaining ranks handle
        the PP work. Otherwise, all ranks do both PP and PME work.

Running mdrun within a single node
----------------------------------

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
    Alternatively, mdrun will honor the appropriate system
    environment variable (e.g. ``OMP_NUM_THREADS``) if set.

``-npme``
    The total number of ranks to dedicate to the long-ranged
    component of PME, if used. The default, -1, will dedicate ranks
    only if the total number of threads is at least 12, and will use
    around a quarter of the ranks for the long-ranged component.

``-ntomp_pme``
    When using PME with separate PME ranks,
    the total number of OpenMP threads per separate PME ranks.
    The default, 0, copies the value from ``-ntomp``.

``-gpu_id``
    A string that specifies the ID numbers of the GPUs to be
    used by corresponding PP ranks on this node. For example,
    "0011" specifies that the lowest two PP ranks use GPU 0,
    and the other two use GPU 1.

``-pin``
    Can be set to "auto," "on" or "off" to control whether
    mdrun will attempt to set the affinity of threads to cores.
    Defaults to "auto," which means that if mdrun detects that all the
    cores on the node are being used for mdrun, then it should behave
    like "on," and attempt to set the affinities (unless they are
    already set by something else).

``-pinoffset``
    If ``-pin on``, specifies the logical core number to
    which mdrun should pin the first thread. When running more than
    one instance of mdrun on a node, use this option to to avoid
    pinning threads from different mdrun instances to the same core.

``-pinstride``
    If ``-pin on``, specifies the stride in logical core
    numbers for the cores to which mdrun should pin its threads. When
    running more than one instance of mdrun on a node, use this option
    to to avoid pinning threads from different mdrun instances to the
    same core.  Use the default, 0, to minimize the number of threads
    per physical core - this lets mdrun manage the hardware-, OS- and
    configuration-specific details of how to map logical cores to
    physical cores.

``-ddorder``
    Can be set to "interleave," "pp_pme" or "cartesian."
    Defaults to "interleave," which means that any separate PME ranks
    will be mapped to MPI ranks in an order like PP, PP, PME, PP, PP,
    PME, ... etc. This generally makes the best use of the available
    hardware. "pp_pme" maps all PP ranks first, then all PME
    ranks. "cartesian" is a special-purpose mapping generally useful
    only on special torus networks with accelerated global
    communication for Cartesian communicators. Has no effect if there
    are no separate PME ranks.

``-nb``
    Used to set where to execute the non-bonded interactions.
    Can be set to "auto", "cpu", "gpu."
    Defaults to "auto," which uses a compatible GPU if available.
    Setting "cpu" requires that no GPU is used. Setting "gpu" requires
    that a compatible GPU be available and will be used.

Examples for mdrun on one node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    gmx mdrun

Starts mdrun using all the available resources. mdrun
will automatically choose a fairly efficient division
into thread-MPI ranks, OpenMP threads and assign work
to compatible GPUs. Details will vary with hardware
and the kind of simulation being run.

::

    gmx mdrun -nt 8

Starts mdrun using 8 threads, which might be thread-MPI
or OpenMP threads depending on hardware and the kind
of simulation being run.

::

    gmx mdrun -ntmpi 2 -ntomp 4

Starts mdrun using eight total threads, with four thread-MPI
ranks and two OpenMP threads per core. You should only use
these options when seeking optimal performance, and
must take care that the ranks you create can have
all of their OpenMP threads run on the same socket.
The number of ranks must be a multiple of the number of
sockets, and the number of cores per node must be
a multiple of the number of threads per rank.

::

    gmx mdrun -gpu_id 12

Starts mdrun using GPUs with IDs 1 and 2 (e.g. because
GPU 0 is dedicated to running a display). This requires
two thread-MPI ranks, and will split the available
CPU cores between them using OpenMP threads.

::

    gmx mdrun -ntmpi 4 -gpu_id "1122"

Starts mdrun using four thread-MPI ranks, and maps them
to GPUs with IDs 1 and 2. The CPU cores available will
be split evenly between the ranks using OpenMP threads.

::

    gmx mdrun -nt 6 -pin on -pinoffset 0
    gmx mdrun -nt 6 -pin on -pinoffset 3

Starts two mdrun processes, each with six total threads.
Threads will have their affinities set to particular
logical cores, beginning from the logical core
with rank 0 or 3, respectively. The above would work
well on an Intel CPU with six physical cores and
hyper-threading enabled. Use this kind of setup only
if restricting mdrun to a subset of cores to share a
node with other processes.

::

    mpirun -np 2 gmx_mpi mdrun

When using an :ref:`gmx mdrun` compiled with external MPI,
this will start two ranks and as many OpenMP threads
as the hardware and MPI setup will permit. If the
MPI setup is restricted to one node, then the resulting
:ref:`gmx mdrun` will be local to that node.

Running mdrun on more than one node
-----------------------------------
This requires configuring |Gromacs| to build with an external MPI
library. By default, this mdrun executable is run with
:ref:`mdrun_mpi`. All of the considerations for running single-node
mdrun still apply, except that ``-ntmpi`` and ``-nt`` cause a fatal
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
    mdrun features are not compatible with this, and these ignore
    this option.

``-dlb``
    Can be set to "auto," "no," or "yes."
    Defaults to "auto." Doing Dynamic Load Balancing between MPI ranks
    is needed to maximize performance. This is particularly important
    for molecular systems with heterogeneous particle or interaction
    density. When a certain threshold for performance loss is
    exceeded, DLB activates and shifts particles between ranks to improve
    performance.

``-gcom``
    During the simulation :ref:`gmx mdrun` must communicate between all ranks to
    compute quantities such as kinetic energy. By default, this
    happens whenever plausible, and is influenced by a lot of [.mdp
    options](#mdp-options). The period between communication phases
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

Examples for mdrun on more than one node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The examples and explanations for for single-node mdrun are
still relevant, but ``-nt`` is no longer the way
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

    mpirun -np 4 gmx mdrun -ntomp 6 -gpu_id 00

Starts :ref:`mdrun_mpi` on a machine with two nodes, using
four total ranks, each rank with six OpenMP threads,
and both ranks on a node sharing GPU with ID 0.

::

    mpirun -np 8 gmx mdrun -ntomp 3 -gpu_id 0000

Using a same/similar hardware as above,
starts :ref:`mdrun_mpi` on a machine with two nodes, using
eight total ranks, each rank with three OpenMP threads,
and all four ranks on a node sharing GPU with ID 0.
This may or may not be faster than the previous setup
on the same hardware.

::

    mpirun -np 20 gmx_mpi mdrun -ntomp 4 -gpu_id 0

Starts :ref:`mdrun_mpi` with 20 ranks, and assigns the CPU cores evenly
across ranks each to one OpenMP thread. This setup is likely to be
suitable when there are ten nodes, each with one GPU, and each node
has two sockets.

::

    mpirun -np 20 gmx_mpi mdrun -gpu_id 00

Starts :ref:`mdrun_mpi` with 20 ranks, and assigns the CPU cores evenly
across ranks each to one OpenMP thread. This setup is likely to be
suitable when there are ten nodes, each with one GPU, and each node
has two sockets.

::

    mpirun -np 20 gmx_mpi mdrun -gpu_id 01

Starts :ref:`mdrun_mpi` with 20 ranks. This setup is likely
to be suitable when there are ten nodes, each with two
GPUs.

::

    mpirun -np 40 gmx_mpi mdrun -gpu_id 0011

Starts :ref:`mdrun_mpi` with 40 ranks. This setup is likely
to be suitable when there are ten nodes, each with two
GPUs, and OpenMP performs poorly on the hardware.

Controlling the domain decomposition algorithm
----------------------------------------------
This section lists all the options that affect how the domain
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

Finding out how to run mdrun better
-----------------------------------

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
* Born radii
* Force
* Waiting + Communication of force
* Particle mesh Ewald
* PME redist. X/F
* PME spread/gather
* PME 3D-FFT
* PME 3D-FFT Communication
* PME solve Lennard-Jones
* PME solve Elec
* PME wait for particle-particle
* Wait + Receive PME force
* Wait GPU nonlocal
* Wait GPU local
* Non-bonded position/force buffer operations
* Virtual site spread
* COM pull force
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

Counters will appear in MD log file only if the related parts of the code were
executed during the :ref:`gmx mdrun` run. There is also a special counter called "Rest" which
indicated for the amount of time not accounted for by any of the counters above. Theerfore,
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
* Nonbonded force
* Ewald force correction
* Non-bonded position buffer operations
* Non-bonded force buffer operations

Subcounters are geared toward developers and have to be enabled during compilation. See
:doc:`/dev-manual/build-system` for more information.

TODO In future patch:
- red flags in log files, how to interpret wallcycle output
- hints to devs how to extend wallcycles

TODO In future patch: import wiki page stuff on performance checklist; maybe here,
maybe elsewhere

Running mdrun with GPUs
-----------------------

NVIDIA GPUs from the professional line (Tesla or Quadro) starting with
the Kepler generation (compute capability 3.5 and later) support changing the
processor and memory clock frequency with the help of the applications clocks feature.
With many workloads, using higher clock rates than the default provides significant
performance improvements.
For more information see the `NVIDIA blog article`_ on this topic.
For |Gromacs| the highest application clock rates are optimal on all hardware
available to date (up to and including Maxwell, compute capability 5.2).

Application clocks can be set using the NVIDIA system managemet tool
``nvidia-smi``. If the system permissions allow, :ref:`gmx mdrun` has
built-in support to set application clocks if built with NVML support. # TODO add ref to relevant section
Note that application clocks are a global setting, hence affect the
performance of all applications that use the respective GPU(s).
For this reason, :ref:`gmx mdrun` sets application clocks at initialization
to the values optimal for |Gromacs| and it restores them before exiting
to the values found at startup, unless it detects that they were altered
during its runtime.

.. _NVIDIA blog article: https://devblogs.nvidia.com/parallelforall/increase-performance-gpu-boost-k80-autoboost/

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
In general, there a user can do little to avoid such overheads, but there
are a few cases where tweaks can give performance benefits.
In single-rank runs timing of GPU tasks is by default enabled and,
while in most cases its impact is small, in fast runs performance can be affected.
The performance impact will be most significant on NVIDIA GPUs with CUDA,
less on AMD with OpenCL.
In these cases, when more than a few percent of "Launch GPU ops" time is observed,
it is recommended turning off timing by setting the ``GMX_DISABLE_GPU_TIMING``
environment variable.
In parallel runs with with many ranks sharing a GPU
launch overheads can also be reduced by staring fewer thread-MPI
or MPI ranks per GPU; e.g. most often one rank per thread or core is not optimal.

The second type of overhead, interference of the GPU driver with CPU computation,
is caused by the scheduling and coordination of GPU tasks.
A separate GPU driver thread can require CPU resources
which may clash with the concurrently running non-offloaded tasks,
potentially degrading the performance of PME or bonded force computation.
This effect is most pronounced when using AMD GPUs with OpenCL with
all stable driver releases to date (up to and including fglrx 12.15).
To minimize the overhead it is recommended to
leave a CPU hardware thread unused when launching :ref:`gmx mdrun`,
especially on CPUs with high core count and/or HyperThreading enabled.
E.g. on a machine with a 4-core CPU and eight threads (via HyperThreading) and an AMD GPU,
try ``gmx mdrun -ntomp 7 -pin on``.
This will leave free CPU resources for the GPU task scheduling
reducing interference with CPU computation.
Note that assigning fewer resources to :ref:`gmx mdrun` CPU computation
involves a tradeoff which may outweigh the benefits of reduced GPU driver overhead,
in particular without HyperThreading and with few CPU cores.

TODO In future patch: any tips not covered above

Running the OpenCL version of mdrun
-----------------------------------

The current version works with GCN-based AMD GPUs, and NVIDIA CUDA
GPUs. Make sure that you have the latest drivers installed. For AMD GPUs,
Mesa version 17.0 or newer with LLVM 4.0 or newer is supported in addition
to the proprietary driver. For NVIDIA GPUs, using the proprietary driver is
required as the open source nouveau driver (available in Mesa) does not
provide the OpenCL support.
The minimum OpenCL version required is |REQUIRED_OPENCL_MIN_VERSION|. See
also the :ref:`known limitations <opencl-known-limitations>`.

Devices from the AMD GCN architectures (all series) and NVIDIA Fermi
and later (compute capability 2.0) are known to work, but before
doing production runs always make sure that the |Gromacs| tests
pass successfully on the hardware.

The OpenCL GPU kernels are compiled at run time. Hence,
building the OpenCL program can take a few seconds introducing a slight
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

- No Intel devices (CPUs, GPUs or Xeon Phi) are supported
- Due to blocking behavior of some asynchronous task enqueuing functions
  in the NVIDIA OpenCL runtime, with the affected driver versions there is
  almost no performance gain when using NVIDIA GPUs.
  The issue affects NVIDIA driver versions up to 349 series, but it
  known to be fixed 352 and later driver releases.
- On NVIDIA GPUs the OpenCL kernels achieve much lower performance
  than the equivalent CUDA kernels due to limitations of the NVIDIA OpenCL
  compiler.
- The AMD APPSDK version 3.0 ships with OpenCL compiler/runtime components,
  libamdocl12cl64.so and libamdocl64.so (only in earlier releases),
  that conflict with newer fglrx GPU drivers which provide the same libraries.
  This conflict manifests in kernel launch failures as, due to the library path
  setup, the OpenCL runtime loads the APPSDK version of the aforementioned
  libraries instead of the ones provided by the driver installer.
  The recommended workaround is to remove or rename the APPSDK versions of the
  offending libraries.

Limitations of interest to |Gromacs| developers:

- The current implementation is not compatible with OpenCL devices that are
  not using warp/wavefronts or for which the warp/wavefront size is not a
  multiple of 32
- Some Ewald tabulated kernels are known to produce incorrect results, so
  (correct) analytical kernels are used instead.

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
* Use a new compiler version, especially for gcc (e.g. from the version 5 to 6
  the performance of the compiled code improved a lot).
* MPI library: OpenMPI usually has good performance and causes little trouble.
* Make sure your compiler supports OpenMP (some versions of Clang don't).
* If you have GPUs that support either CUDA or OpenCL, use them.

  * Configure with ``-DGMX_GPU=ON`` (add ``-DGMX_USE_OPENCL=ON`` for OpenCL).
  * For CUDA, use the newest CUDA availabe for your GPU to take advantage of the
    latest performance enhancements.
  * Use a recent GPU driver.
  * If compiling on a cluster head node, make sure that ``GMX_CPU_ACCELERATION``
    is appropriate for the compute nodes.

Run setup
^^^^^^^^^

* For an approximately spherical solute, use a rhombic dodecahedron unit cell.
* When using a time-step of 2 fs, use :mdp:`cutoff-scheme` = :mdp:`h-bonds`
  (and not :mdp:`all-bonds`), since this is faster, especially with GPUs,
  and most force fields have been parametrized with only bonds involving
  hydrogens constrained.
* You can increase the time-step to 4 or 5 fs when using virtual interaction
  sites (``gmx pdb2gmx -vsite h``).
* For massively parallel runs with PME, you might need to try different numbers
  of PME ranks (``gmx mdrun -npme ???``) to achieve best performance;
  ``gmx tune_pme`` can help automate this search.
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
