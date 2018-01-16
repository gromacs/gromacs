Bugs fixed
^^^^^^^^^^

Fixed multiple time stepping with Parrinello-Rahman and Nose-Hoover.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
These now work in correct Trotter style, applied once and scaled by
the correct number of steps.

Fixes :issue:`2031`
Fixes :issue:`2032`

Applied Berendsen pressure coupling only at nstpcouple steps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Berendsen pressure coupling was mistakenly applied on successive
steps. Since there is no need for this, this is changed to act only on
nstpcouple steps. Note that this change prevents continuation from old
checkpoint files for Berendsen pressuring-coupling runs, since the
previous-step pressure is no longer stored.

Added missing Ewald correction for PME-User
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With :mdp-value:`coulombtype=PME-User`, the Ewald mesh energy was not subtracted
leading to (very) incorrect Coulomb energies and forces.

Fixes :issue:`2286`

Fixed incorrect dV/dlambda for walls
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The free-energy derivative dV/dlambda for walls, which can be
perturbed by changing atom types of non-wall atoms, only contained the
B-state contribution.

Fixes :issue:`2267`

Supported OpenMP for orientation restraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Previously this was broken, but has been fixed and is now tested and
supported.

Fixed orientation restraint reference
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The resetting of the COM of the molecule with orientation restraints
for fitting to the reference structure was done with the COM of the
reference structure instead of the instantaneous structure. This does
not affect the restraining (unless ensemble averaging is used), only
the printed orientation tensor.

Fixes :issue:`2219`

Used graph with orientation restraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With the Verlet cut-off scheme by default molecules are not made whole.
Now they are made whole when orientation restraints are used.
Added checks and assertions for correct PBC treatment with orientation
restraints.

Fixes :issue:`2228`

Fixed Ekin at step 0 with COM removal
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The kinetic energy at step 0 was computed from the velocities without
the center of mass velocity removed. This could cause a relatively
large jump in kinetic energy, especially for small systems.
Now compute_globals is called twice with COM removal so we get
the correct kinetic energy.

Fixed :ref:`gmx grompp` with Andersen massive and no COM removal
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixed a floating point exception leading to a segv.
Also fixed possible different rounding for the interval for
Andersen massive in :ref:`gmx grompp` in mdrun for the common case where tau-t
is a multiple of delta-t.

Fixes :issue:`2256`

Improved Verlet buffer constraint estimate
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The displacement estimate for a constrained atom (typically H)
rotating around the COM with a partner atom is now derived and
documented correctly.  Note that we (still) use a Gaussian with
matched variance, which results in a much larger buffer than
necessary, since the tail of the displacement distribution sets the
buffer size and the Gaussian has a long tail whereas the actual
distribution has no tail.

Fixed virtual site generation for water oxygens not named OW
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:ref:`gmx pdb2gmx` would break when generating virtual sites if water oxygens
were not named OW. Now checking for the atomnumber instead.

Fixes :issue:`2268`

Fixed thread-MPI rank choice for orientation restraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Only a single rank is supported, so that must be what the thread-MPI
code will choose. There's another check later on that catches the
multi-rank MPI case.

Fixed some incorrect behavior with :ref:`gmx solvate`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:ref:`gmx solvate` cannot replicate non-rectangular solvent boxes correctly
(there are several different places that assume a diagonal box matrix),
so give a fatal error if that is attempted.  To support some uses with
triclinic boxes, skip the replication step if the solvent and target box
sizes are already equal.

Support for general triclinic boxes can be added separately, and the
check introduced here can be valuable even in that case: it keeps a
pre-equilibrated solvent box intact if the target box size is the same.

Related to fix of :issue:`2148`

Fixed DD exact continuation in reproducible node
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With domain decomposition, the local atom density, used for setting
the search grid for sorting particles, was based on the local atom
count including atoms/charge groups that would be moved to neighboring
cells. This lead to a different density value, which in turn could
result in a different number of search grid cells and thus a different
summation order during a run compared with continuing that run from a
checkpoint, when no atoms would be moved. That difference violated
the intention of ``mdrun -reprod``, and is now fixed.

Refs Fixes :issue:`2318`

Made mdrun only stop at nstlist steps with mdrun -reprod
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Stopping mdrun with two INT or TERM signals (e.g. from Ctrl-C from the
terminal shell) would always happen right after the first global
communication step. But this breaks exact continuation. Now with
``mdrun -reprod`` a second signal will still stop at a pair-list
generation step, like with the first signal, so we can still have
exact continuation.

Fixes :issue:`2318`

Added check for GPU detection support before detecting GPU devices
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
When a CUDA-enabled binary was run on a node with no CUDA driver
available, a note was issued that the version of the CUDA driver is
insufficient, which was wrong and now fixed.

Fixes :issue:`2322`

Removed duplicated lines from OPLS ffbonded.itp
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Identical lines have been removed, as identified with uniq.

Fixes :issue:`1678`.

mdrun no longer warns about NVML clocks that are at max
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If the clocks are already maxed out there is no point in echoing
warnings about not being able to set them.

Fixes :issue:`2313`.

Used reduced default tolerances for tpx comparison
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The tolerances for gmx check are mainly intended for handling slight
statistical deviations, but they can hide differences between tpr
files, when the user likely wants exact checks on small quantities
like Lennard-Jones parameters. This changes changes the default
relative tolerance to 0.000001 and the absolute tolerance to zero, so
that we only allow for any minor differences due to compiler
optimization.

Fixes :issue:`2024`.

Fixed return values of frame-reading functions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This function was based on read_first_x that returned the number of
atoms, and was documented to do the same, but has always returned a
logical boolean about whether a frame has been read. This led to
aspects of ``gmx spatial`` and ``gmx trjcat -demux`` being broken.

Fixed by returning a proper bool, and fixing the remaining logic that
used the return value in a non-boolean sense.

Refs :issue:`2157`

Removed PBC before generating TPR with group scheme
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Ensure that all molecules have been made whole before generating the
run input file when using the group scheme, to avoid error messages
for large charge groups when molecules are broken over PBC boundaries.

Fixes :issue:`2339`

Fixed PBC error in gmx_spatial
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Fixes :issue:`2157`.

Documented power spectrum options of gmx velacc
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2019`.

Changed to require .tpr file for gmx cluster
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The program could crash without it, so it wasn't optional.

Fixes :issue:`2170`.

Disallowed ascii formats for gmx trjcat
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Since gmx trjcat (deliberately) does not use any .tpr file, the tool
can't handle trajectory formats such as .gro or .pdb where
atom/residue names are needed.

Fixes :issue:`2225`.

Improved grompp missing-parameters error message
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If an interaction entry had parameters but not the function type, then
the error message has been confusing. Note that even when only one
function type is implemented, the field is still required, which makes
for ready extensibility.

Refs :issue:`2144`

Checked for large energy at first step
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Also added step number to fatal error message.

Fixes :issue:`2333`

Disallowed combination of PME-user and verlet cutoff
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2332`

Avoided confusing message at end of non-dynamical runs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Energy minimization, test-particle insertion, normal-mode analysis,
etc.  are not targets for performance optimization so we will not
write performance reports. This commit fixes an oversight whereby we
would warn a user when the lack of performance report is normal and
expected.

Fixes :issue:`2172`

Changed to require ``-ntmpi`` when setting ``-ntomp`` and using GPUs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With GPUs and thread-MPI, setting only ``gmx mdrun -ntomp`` could lead
to oversubscription of the hardware threads.  Now, with GPUs and
thread-MPI the user is required to set ``-ntmpi`` when using
``-ntomp``. Here we chose that to also require ``-ntmpi`` when the
user specified both ``-nt`` and ``-ntomp``; here we could infer the
number of ranks, but it's safer to ask the user to explicity set
``-ntmpi``.  Note that specifying both ``-ntmpi`` and ``-nt`` has
always worked correctly.

Fixes :issue:`2348`
