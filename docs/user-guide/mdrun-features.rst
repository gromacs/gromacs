Useful mdrun features
=======================
This section discusses features in :ref:`gmx mdrun` that don't fit well
elsewhere.

.. _single-point energy:

Re-running a simulation
-----------------------
The rerun feature allows you to take any trajectory file ``traj.trr``
and compute quantities based upon the coordinates in that file using
the model physics supplied in the ``topol.tpr`` file. It can be used
with command lines like ``mdrun -s topol -rerun traj.trr``. That :ref:`tpr`
could be different from the one that generated the trajectory. This
can be used to compute the energy or forces for exactly the
coordinates supplied as input, or to extract quantities based on
subsets of the molecular system (see :ref:`gmx convert-tpr` and
:ref:`gmx trjconv`). It is easier to do a correct "single-point" energy
evaluation with this feature than a 0-step simulation.

Neighbour searching is normally performed for every frame in the
trajectory, since :ref:`gmx mdrun` can no longer assume anything about how the
structures were generated. If :mdp:`nstlist` is zero, then only one
neighbour list will be constructed. Naturally, no update or constraint
algorithms are ever used.

Running a simulation in reproducible mode
-----------------------------------------
It is generally difficult to run an efficient parallel MD simulation
that is based primarily on floating-point arithmetic and is fully
reproducible. By default, :ref:`gmx mdrun` will observe how things are going
and vary how the simulation is conducted in order to optimize
throughput. However, there is a "reproducible mode" available with
``mdrun -reprod`` that will systematically eliminate all sources of
variation within that run; repeated invocations on the same input and
hardware will be binary identical. However, running in this mode on
different hardware, or with a different compiler, etc. will not be
reproducible. This should normally only be used when investigating
possible problems.

Running multi-simulations
-------------------------
There are numerous situations where running a related set of
simulations within the same invocation of mdrun are necessary or
useful. Running a replica-exchange simulation requires it, as do
simulations using ensemble-based distance or orientation restraints.
Running a related series of lambda points for a free-energy
computation is also convenient to do this way.

This feature requires
`configuring GROMACS with an external MPI library <../install-guide/main.html#mpi-support>`_
so that the set of
simulations can communicate. The ``n`` simulations within the set can
use internal MPI parallelism also, so that ``mpirun -np x mdrun_mpi``
for ``x`` a multiple of ``n`` will use ``x/n`` ranks per simulation.

There are two ways of organizing files when running such
simulations. All of the normal mechanisms work in either case,
including ``-deffnm``.

``-multidir``
   You must create a set of ``n`` directories for the ``n`` simulations,
   place all the relevant input files in those directories (e.g. named
   ``topol.tpr``), and run with
   ``mpirun -np x mdrun_mpi -s topol -multidir <names-of-directories>``.
   If the order of the simulations
   within the multi-simulation is significant, then you are responsible
   for ordering their names when you provide them to ``-multidir``. Be
   careful with shells that do filename globbing dictionary-style, e.g.
   ``dir1 dir10 dir11 ... dir2 ...``. This option is generally the
   most convenient to use. ``mdrun -table`` for the group cutoff-scheme
   works only in this mode.

``-multi``
   You must organize that the filenames for each simulation in a set of
   ``n`` simulations have an integer ``0`` through ``n-1`` appended to
   the filename (e.g. ``topol2.tpr``), and run with
   ``mpirun -np x mdrun -multi n -s input``. The order of simulations
   within the set is determined by the integer appended.

Examples running multi-simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    mpirun -np 32 mdrun_mpi -multi

Starts a multi-simulation on 32 ranks with as many simulations ``n`` as
there are files named ``topol*.tpr`` for integers ``0`` to ``n-1``. Other
input and output files are suffixed similarly.

::

    mpirun -np 32 mdrun_mpi -multidir a b c d

Starts a multi-simulation on 32 ranks with 4 simulations. The input
and output files are found in directories ``a``, ``b``, ``c``, and ``d``.

::

    mpirun -np 32 mdrun_mpi -multidir a b c d -gpu_id 0000000011111111

Starts the same multi-simulation as before. On a machine with two
physical nodes and two GPUs per node, there will be 16 MPI ranks per
node, and 8 MPI ranks per simulation. The 16 MPI ranks doing PP work
on a node are mapped to the GPUs with IDs 0 and 1, even though they
come from more than one simulation. They are mapped in the order
indicated, so that the PP ranks from each simulation use a single
GPU. However, the order ``0101010101010101`` could run faster.

Running replica-exchange simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running a multi-simulation, using ``mdrun -replex n`` means that a
replica exchange is attempted every given number of steps. The number
of replicas is set with the ``-multi`` or ``-multidir`` option, described
above.  All run input files should use a different value for the
coupling parameter (e.g. temperature), which ascends over the set of
input files. The random seed for replica exchange is set with
``-reseed``. After every exchange, the velocities are scaled and
neighbor searching is performed. See the Reference Manual for more
details on how replica exchange functions in GROMACS.
