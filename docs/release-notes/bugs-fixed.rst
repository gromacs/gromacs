.. TODO Remove beta-phase fixes below before final release

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

Add missing Ewald correction for PME-User
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With :mdp-value:`coulombtype=PME-User`, the Ewald mesh energy was not subtracted
leading to (very) incorrect Coulomb energies and forces.

Fixes :issue:`2286`

Fix incorrect dV/dlambda for walls
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The free-energy derivative dV/dlambda for walls, which can
be perturbed by changing atom types of non-wall atoms, only
contained the B-state contribution.

Fixes :issue:`2267`

Supported OpenMP for orientation restraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Previously this was broken, but has been fixed and is now tested
and supported.

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

Fix Ekin at step 0 with COM removal
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

Now mdrun only stops at nstlist steps with mdrun -reprod
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
Identical lines have been removed, as identified
with uniq.

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

``mdrun -pme cpu -pmefft gpu`` now gives a fatal error  - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Previously was silently ignored.

Fixed mdrun -nb auto -pme auto when GPUs are absent - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The logic was flawed such that GPUs were "selected" for use even
though none had been detected. That led to the GPU behaviour of
avoiding using separate PME ranks.

Also made a minor fix to the logic for emulation. The new
interpretation of ``mdrun -gpu_id`` does not need to trigger an error
when GPU IDs have been supplied along with the emulation environmnet
variable.

Fixes :issue:`2315`

Fixed ArrayRef<SimdDInt32> for SSE/AVX128 - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2326`

Fixed PME gather in double with AVX(2)_128 - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The 4NSIMD PME gather change did not change the conditional
for grid alignment. This is made consistent here.
Note that the 4NSIMD change lowered the performance of PME gather
on AVX_128_FMA and AVX2_128 in double precision. We should consider
using 256-bit AVX for double precision instead.

Fixes :issue:`2326`

Reformulated PME and SHAKE test tolerances - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2306`
Fixes :issue:`2337`
Fixes :issue:`2338`

Fixed freeing of GPU context - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If a device context was not used, CUDA gives an error if we attempt to
clear it, so we must avoid clearing it.

:issue:`2322`

Fixed initial temperature reporting - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixes :issue:`2314`

Fixed compilation issues for AVX-512 - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
- gcc-5.4.0 incorrectly requires the second argument of
  _mm512_i32gather_pd() to be a double pointer instead
  of void, but this should fix compilation for both
  cases.
- Work around double precision permute instruction
  only available with AVX512VL instructions.

Fixes :issue:`2312`

Cleared vsite velocities for simple integrators - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The simple integrator loops did not clear
the velocities of virtual sites. This allows velocities of virtual
sites to slowly increase over time. To prevent this, velocities
of virtual sites are now cleared in a separate loop.

Fixes :issue:`2316`

Fixed fft5d pinning - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A CUDA build on a node with no driver installed can never have
selected a CUDA pinning policy, and erroneously unpinning leads to a
fatal error. Instead, FFT5D now remembers whether it made pinning
possible, which can only occur when there was a driver and a valid
device, so that it can unpin only when appropriate.

Fixes :issue:`2322`

Avoided assertion failure in AWH - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
With an unstable reaction coordinate or unequilibrated system, AWH
could cause an assertion to fail. Now AWH checks for valid coordinate
input and throws an exception with a clear message.

Corrected AWH input file name in documentation - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Mdrun was expecting user input data file 'awhinit.xvg' while the
mdp-option documentation has 'awh-init.xvg'.

Changed the GPU SMT cut-off to quadratic - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The advantage of SMT diminishes rapidly with the number of cores.
So the system sizes should be compares to the square of the number
of cores.

Fixed AVX-512 SIMD test for C compilation - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Avoid using C++ features in the test, since it should test both the C
and C++ compilers.

Leave NVML use off by default - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Even if NVML is found, leave the default off because the
linking is unreliable for reasons that are currently unclear,
and only in some cases is linking with NVML advantageous.

Fixes :issue:`2311`

Fixes for compiler support - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Double precision, debug mode, proper release mode and some quirky
cases were all improved in multiple ways to compile and pass tests
reliably.

Consume any error produced during GPU detection - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Having reported it, we should clear the CUDA error status so that
future calls do not continue to return it.

Fixes :issue:`2321`

Replace intrinsic with inline asm for AVX512 unit test - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Using inline assembly avoids compilers at low optimization
levels not generating efficient code for the timing routines, and
also avoids needing an assembler.

Fixes :issue:`2340`

Fixed table tests and improve table construction - beta phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Since compilers are allowed to use different FMA constructs, we
now allow the consistency check to deviate a few ulps.

For sinc and other extreme functions that oscillate, the
scan over the definition range to locate the minimum quotient
between the 1st and 4th derivative to set the table spacing
exposes some delicate errors. Basically, it is not possible
to have arbitrarily low relative errors for the derivative
for a function that has large magnitude in the same place.
For now we reduce the test interval for sinc(); this should
anyway not be relevant for normal well-behaved MD functional
forms.

Fixes :issue:`2336`.

Supported Simd4N for SimdRealWidth<4 - beta-phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If the SIMD with is smaller 4 but Simd4N is supported
then use Simd4 for Simd4N.

Fixes :issue:`2327`

Made AVX-512 CMake detection work - beta phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Both inline assembly and the support flag have to be set for the
timing code to be compiled.

Fixed shift usage for KNC - beta phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
9437181eacb removed the shift operator without replacing the usage for
KNC.

Made acceleration correction VCM mode work - beta phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The new acceleration correction VCM mode did not actually correct
the coordinate for the acceleration, since a null pointer was passed.
Introduced an extra CGLO flag to allow for correction of the
coordinates, but leave the initial coordinates unaffected.

Fix builds on ARM & clarify (ARM) GPU support - beta phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixed a typo in architecture.h that prevented
the Neon Asimd instructions from being selected,
and updated the CPU brand detection to also look
for a new label with Tegra X1 on Ubuntu 16.04

Fixes :issue:`2287`

Improved documentation and code for physical validation - beta phase fix
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Adds documentation for the physical validation suite in
docs/dev-manual/physical_validation.rst

As this was misunderstandable, changed the default behavior of
`make check-phys` and `make check-all` to actually run the simulations.
This might take very long, but since the physical validation tests need to
be turned on explicitly via cmake option, the chances of somebody using the
tests by mistake are low. The `check` targets are:

* `make check`: Run unit and regression tests (unchanged)
* `make check-phys`: Run simulations needed for physical validation, then
  run physical validation tests
* `make check-phys-analyze`: Only run physical validation tests, assuming
  that simulations were run previously and are available.
* `make check-all`: Combination of `make check` and `make check-phys`

Additionally, `make check-phys-prepare` can be used to prepare |Gromacs|
input files and a script to run the simulations needed for the physical
validation tests.

Fixes :issue:`2349`
