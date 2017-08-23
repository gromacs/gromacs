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

Neighbor searching is normally performed for every frame in the
trajectory, since :ref:`gmx mdrun` can no longer assume anything about how the
structures were generated. If :mdp:`nstlist` is zero, then only one
neighbor list will be constructed. Naturally, no update or constraint
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
:ref:`configuring |Gromacs| with an external MPI library <mpi-support>`
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
   ``mpirun -np x gmx_mpi mdrun -s topol -multidir <names-of-directories>``.
   If the order of the simulations
   within the multi-simulation is significant, then you are responsible
   for ordering their names when you provide them to ``-multidir``. Be
   careful with shells that do filename globbing dictionary-style, e.g.
   ``dir1 dir10 dir11 ... dir2 ...``. This option is generally the
   most convenient to use. ``gmx mdrun -table`` for the group cutoff-scheme
   works only in this mode.

``-multi``
   You must organize that the filenames for each simulation in a set of
   ``n`` simulations have an integer ``0`` through ``n-1`` appended to
   the filename (e.g. ``topol2.tpr``), and run with
   ``mpirun -np x gmx mdrun -multi n -s input``. The order of simulations
   within the set is determined by the integer appended.

Examples running multi-simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    mpirun -np 32 gmx_mpi mdrun -multi

Starts a multi-simulation on 32 ranks with as many simulations ``n`` as
there are files named ``topol*.tpr`` for integers ``0`` to ``n-1``. Other
input and output files are suffixed similarly.

::

    mpirun -np 32 gmx_mpi mdrun -multidir a b c d

Starts a multi-simulation on 32 ranks with 4 simulations. The input
and output files are found in directories ``a``, ``b``, ``c``, and ``d``.

::

    mpirun -np 32 gmx_mpi mdrun -multidir a b c d -gputasks 0000000011111111

Starts the same multi-simulation as before. On a machine with two
physical nodes and two GPUs per node, there will be 16 MPI ranks per
node, and 8 MPI ranks per simulation. The 16 MPI ranks doing PP work
on a node are mapped to the GPUs with IDs 0 and 1, even though they
come from more than one simulation. They are mapped in the order
indicated, so that the PP ranks from each simulation use a single
GPU. However, the order ``0101010101010101`` could run faster.

Running replica-exchange simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running a multi-simulation, using ``gmx mdrun -replex n`` means that a
replica exchange is attempted every given number of steps. The number
of replicas is set with the ``-multi`` or ``-multidir`` option, described
above.  All run input files should use a different value for the
coupling parameter (e.g. temperature), which ascends over the set of
input files. The random seed for replica exchange is set with
``-reseed``. After every exchange, the velocities are scaled and
neighbor searching is performed. See the Reference Manual for more
details on how replica exchange functions in |Gromacs|.

Controlling the length of the simulation
----------------------------------------

Normally, the length of an MD simulation is best managed through the
[.mdp] option [nsteps](#nsteps), however there are situations where
more control is useful. `gmx mdrun -nsteps 100` overrides the [.mdp] file
and executes 100 steps. `gmx mdrun -maxh 2.5` will terminate the
simulation shortly before 2.5 hours elapse, which can be useful when
running under cluster queues (as long as the queuing system does not
ever suspend the simulation).

Running a membrane protein embedding simulation
-----------------------------------------------

This is a module to help embed a membrane protein into an equilibrated
lipid bilayer at a position and orientation specified by the user. 

This method was initially described as a ProtSqueeze technique 
(`Yesylevskyy S.O., J Chem Inf Model 47(5) (2007) 1986-94`_) and 
later implemented in |Gromacs| as g_membed tool (`Wolf et al, J Comp Chem 31 (2010) 2169-2174`_). 
Currently the functionality of g_membed is available in mdrun if 
``-membed`` option is specified (see below).

.. _Yesylevskyy S.O., J Chem Inf Model 47(5) (2007) 1986-94: https://dx.doi.org/10.1021/ci600553y
.. _Wolf et al, J Comp Chem 31 (2010) 2169-2174: http://onlinelibrary.wiley.com/doi/10.1002/jcc.21507/full

The main advantage is that it is possible to use very complex lipid bilayers
with a number of different components that have been relaxed for a
long time in a previous simulation. In theory that could be accomplished
with a procedure similar to :ref:`gmx solvate`, but since lipids are much larger
than water molecules it will lead to a large vacuum layer between the
protein and membrane if we remove all molecules where any atom is
overlapping. Instead, this module works by first artificially shrinking
the protein in the xy-plane, then it removes lipids that overlap with
a much smaller core, after which we gradually push the protein atoms
back to their initial positions, while using normal dynamics for the
rest of the system so lipids adapt to the protein.

To use membrane embedding, start by building a lipid bilayer that is
just-so-slightly larger in the xy-plane than what you expect to need
in the end, and make sure you have enough water outside the membrane
to accommodate globular domains. Place the protein in the same coordinate
file (and topology) as your lipid bilayer, and make sure it is in the
orientation and position you want right in the middle of the bilayer.

The first settings have to be entered in the mdp file that controls
your simulation. You need an energy group corresponding to your
protein, this group should be frozen (all dimensions), and we should
exclude all interactions inside the protein to avoid problems when it
is distorted. For instance:

::

    integrator     = md
    energygrps     = Protein
    freezegrps     = Protein
    freezedim      = Y Y Y
    energygrp_excl = Protein Protein

You will also need a number of settings for the actual membrane
embedding process. These are entered as similar name and value pairs,
but in the separate text data file ``embed.dat`` that you provide as
the argument to the ``-membed`` option (we refer to the below
when explaining the process). The embedding works in for stages:

1. The protein is resized around its center of mass by a factor
   ``xy`` in the xy-plane (the bilayer plane), and a factor ``z``
   along the z-axis (normal to the bilayer). If the height of the
   protein is the same or smaller than the thickness of the
   membrane, a z-fraction larger than 1.0 can prevent the protein
   from being enveloped by the lipids.

2. All lipid and solvent molecules overlapping with the resized
   protein are removed. All interactions inside the protein are
   turned off to prevent numerical issues for small values of the
   scaling fractions.

3. A single md step is performed, where atoms in the rest of the
   system are moved.

4. The resize factors are adjusted by the small amounts
   (1-xy)/nxy and (1-z)/nz, where ``nxy`` and ``nz`` are the
   number of iterations to use.  The resize factor for the xy-plane
   is adjusted first. The resize factor for the z-direction is not
   changed until the xy factor is 1.0 (after ``nxy`` iterations).

5. Steps 3 and 4 are repeated until the protein has again reached
   its original size, i.e. after nxy+nz iterations. After the
   embedding you might still want to perform a short relaxation.

Parameters that can be specified in ``embed.dat``, with default
values that will be used if the setting is omitted:

- ``xyinit`` (0.5) Resize factor for the protein in the xy
  dimension before starting embedding.

- ``xyend`` (1.0) Final resize factor in the xy dimension.

- ``zinit`` (1.0) Resize factor for the protein in the z
   dimension before starting embedding.

- ``zend`` (1.0) Final resize faction in the z dimension.

- ``nxy`` (1000) Number of iteration for the xy dimension.

- ``nz`` (0) Number of iterations for the z dimension.

- ``rad`` (0.22) Probe radius to check for overlap between
  the group to embed and the membrane.

- ``pieces`` (1) Perform piecewise resize. Select parts of the group
  to insert and resize these with respect to their own geometrical center.

- ``asymmetry`` (no) Allow asymmetric insertion, i.e. the number of
  lipids removed from the upper and lower leaflet will not be checked.

- ``ndiff`` (0) Number of lipids that will additionally be removed
  from the lower (negative number) or upper (positive number)
  membrane leaflet.

- ``maxwarn`` (0) Largest number of membed warnings allowed.
