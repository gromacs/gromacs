Bugs fixed
^^^^^^^^^^

These document fixes for issues that have been fixed for the 2016
release, but which have not been back-ported to other branches.

Fixed two problems related to restarts for velocity-Verlet
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The first problem is more serious; in addition to causing problems
with restarts in most cases for velocity-Verlet integrators plus either
Berendsen or v-rescale temperature-coupling algorithms, the
temperature coupling code was called twice. This made the distribution of
kinetic energies too broad (but with the correct average).
Other algorithm combinations were unaffected.

In the second problem, the initial step after restarts with velocity-Verlet
integrators and either Berendsen or v-rescale temperature-coupling algorithms
had too high a pressure because they used an empty virial matrix that
was only filled with MTTK pressure control. The effects of this bug were
very small; it only affected the volume integration for one step on restarts.

:issue:`1883`

Fixed Verlet buffer calculation with nstlist=1
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Under rare circumstances the Verlet buffer calculation code was
called with nstlist=1, which caused a division by zero. The division
by zero is now avoided.
Furthermore, grompp now also determines and prints the Verlet buffer
sizes with nstlist=1, which provider the user information and adds
consistency checks.

:issue:`1993`

Fixed large file issue on 32-bit platforms
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
At some point gcc started to issue a warning instead of a fatal error
for the checking code; fixed to really generate an error now.

:issue:`1834`

Avoided using abort() for fatal errors
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This avoids situations that produce useless core dumps.

:issue:`1866`

Fixed possible division by zero in polarization code
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Avoided numerical overflow with overlapping atoms in Verlet scheme
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The Verlet-scheme kernels did not allow overlapping atoms, even if
they were not interacting (in contrast to the group kernels). Fixed by
clamping the interaction distance so it can not become smaller than
~6e-4 in single and ~1e-18 in double, and when this number is later
multiplied by zero parameters it will not influence forces. The
clamping should never affect normal interactions; mdrun would
previously crash for distances that were this small.  On Haswell, RF
and PME kernels get 3% and 1% slower, respectively.  On CUDA, RF and
PME kernels get 1% and 2% faster, respectively.

:issue:`1958`

Relax pull PBC check
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The check in the pull code for COM distances close to half the box
was too strict for directional pulling. Now dimensions orthogonal
to the pull vector are no longer checked. (The check was actually
not strict enough for directional pulling along x or y in triclinic
units cells, but that is a corner case.)
Furthermore, the direction-periodic hint is now only printed with
geometry direction.

:issue:`1962`

Add detection for ARMv7 cycle counter support
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
ARMv7 requires special kernel settings to allow cycle
counters to be read. This change adds a cmake setting
to enable/disable counters. On all architectures but ARMv7
it is enabled by default, and on ARMv7 we run a small test
program to see if the can be executed successfully. When
cross-compiling to ARMv7 counters will be disabled, but
either choice can be overridden by setting a value for
GMX_CYCLECOUNTERS in cmake.

:issue:`1933`

Introduced fatal error for too few frames in ``gmx dos``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
To prevent ``gmx dos`` from crashing with an incomprehensible error
message when there are too few frames, test for this.

Part of :issue:`1813`

Properly reset CUDA application clocks
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We now store the application clock values we read when starting mdrun
and reset to these values, but only when clocks have not been changed
(by another process) in the meantime.

:issue:`1846`

Fixed replica-exchange debug output to all go to the debug file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
When ``mdrun -debug`` was selected with replica exchange, some of the
order description was printed to mdrun's log file, but it looks like the
actual numbers were being printed to the debug log. This puts them
both in the debug log.

Fixed gmx mdrun -membed to always run on a single rank
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This used to give a fatal error if default thread-MPI mdrun had chosen
more than one rank, but it will now correctly choose to use a single rank.

Fixed issues with using int for number of simulation steps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Mostly we use a 64-bit integer, but we messed up a few
things.

During mdrun -rerun, edr writing complained about the negative step
number, implied it might be working around it, and threatened to
crash, which it can't do. Silenced the complaint during writing,
and reduced the scope of the message when reading.

Fixed TNG wrapper routines to pass a 64-bit integer like they should.

Made various infrastructure use gmx_int64_t for consistency, and noted
where in a few places the practical range of the value stored in such
a type is likely to be smaller. We can't extend the definition of XTC
or TRR, so there is no proper solution available. TNG is already good,
though.

:issue:`2006`

Fixed trr magic-number reading
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The trr header-reading routine returned an "OK" value even if the
magic number was wrong, which might lead to chaotic results
everywhere.  This led to problems if other code (e.g. cpptraj)
mistakenly wrote a wrong-endian trr file, which was then used with
|Gromacs|. (This should never be a thing for XDR files, which are
defined to be big endian, but such code has existed.)

:issue:`1926`

Changed to use only legal characters in OpenCL cache filename
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The option to cache JIT-compiled OpenCL short-ranged kernels needed to
be hardened, so that mdrun would write files whose names would usually
be specific to the device, but also only contain filenames that would
work everywhere, ie only alphanumeric characters from the current
locale.

Fixes for bugs introduced during development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These document fixes for issues that were identified as having been
introduced into the release-2016 branch since it diverged from
release-5-1. These will not appear in the final release notes, because
no formal release is thought to have had the problem. Of course, the
tracked `issues <https://gitlab.com/gromacs/gromacs/-/issues/>`_
remain available should further discussion arise.

Fixed bug in v-rescale thermostat & replica exchange
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Commit 2d0247f6 made random numbers for the v-rescale thermostat that
did not vary over MD steps, and similarly the replica-exchange random
number generator was being reset in the wrong place.

:issue:`1968`

Fixed vsite bug with MPI+OpenMP
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The recent commit b7e4f30d caused non-local virtual sites not be
treated when using OpenMP. This means their coordinates lagged one
step behind and their forces are not spread to the atoms, leading
to small errors in the forces. Note that non-local virtual sites are
only used when local virtual sites use them as a constructing atom;
the most common case is a C/N in a CH3/NH3 group with vsite H's.
Also added a check on the vsite count for debug builds.

:issue:`1981`

Fixed some thread affinity cases
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixed one deadlock in newly refactored thread-affinity code, which
happened with automatic pinning, if only part of the nodes were full.

There is one deadlock still theoretically possible: if thread-MPI
reports that setting the affinity is not possible only on a subset of
ranks, the code deadlocks.  This has always been there and might never
happen, so it is not fixed here.

Removed OpenMP overhead at high parallelization
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Commit 6d98622d introduced OpenMP parallelization for for loops
clearing rvecs of increasing rvecs. For small numbers of atoms per
MPI rank this can increase the cost of the loop by up to a factor 10.
This change disables OpenMP parallelization at low atom count.

Removed std::thread::hardware_concurrency()
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
We should not use std::thread::hardware_concurrency() for determining
the logical processor count, since it only provides a hint.
Note that we still have 3 different sources for this count left.

Added support for linking against external TinyXML-2
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This permits convenient packaging of |Gromacs| by distributions, but
it got lost from gerrit while rebasing.

:issue:`1956`

Fixed data race in hwinfo with thread-MPI
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`1983`

Fixes for Power7 big-endian
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Now compiles and passes all tests in both double and single precision
with gcc 4.9.3, 5.4.0 and 6.1.0 for big-endian VSX.

The change for the code in incrStoreU and decrStoreU addresses an
apparent regression in 6.1.0, where the compiler thinks the type
returned by vec_extract is a pointer-to-float, but my attempts a
reduced test case haven't reproduced the issue.

Added some test cases that might hit more endianness cases in future.

We have not been able to test this on little-endian Power8; there is
a risk the gcc-specific permutations could be endian-sensitive. We'll
test this when we have hardware access, or if somebody runs the tests
for us.

:issue:`1997`
:issue:`1988`

Reduce hwloc & cpuid test requirements
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
On some non-x86 linux platforms hwloc does not report
caches, which means it will fail our strict test
requirements of full topology support. There is no
problem whatsoever with this, so we reduce the
test to only require basic support from hwloc - this
is still better than anything we can get ourselves.
Similarly for CPUID, it is not an error for an
architecture to not provide any of the specific flags
we have defined, so avoid marking it as such.

:issue:`1987`

Work around compilation issue with random test on 32-bit machines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
gcc 4.8.4 running on 32-bit Linux fails a few
tests for random distributions. This seems
to be caused by the compiler doing something
strange (that can lead to differences in the lsb)
when we do not use the result as floating-point
values, but rather do exact binary comparisions.
This is valid C++, and bad behaviour of the
compiler (IMHO), but technically it is not required
to produce bitwise identical results at high
optimization. However, by using floating-point
tests with zero ULP tolerance the problem
appears to go away.

:issue:`1986`

Updated ``gmx wham`` for the new pull setup
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This bring ``gmx wham`` up to date with the new pull setup where the pull
type and geometry can now be set per coordinate and the pull
coordinate has changed and is more configurable.

Fix membed with partial revert of 29943f
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The membrane embedding algorithm must be initialized before
we call init_forcerec(), so it cannot trivially be moved into
do_md(). This has to be cleaned up anyway for release-2017
since we will remove the group scheme be then, but for now
this fix will allow us have the method working in release-2016.

:issue:`1998`
