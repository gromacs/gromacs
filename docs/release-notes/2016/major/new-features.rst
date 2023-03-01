New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

Changed to require a C++11 compiler
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
|Gromacs| now requires both a C++11 and C99 compiler. For details, see
the install guide.

Changed to support only CUDA 5.0 and more recent versions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:issue:`1831`

Allowed rcoulomb > rvdw with PME
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
|Gromacs| has had kernels that support Coulomb PME + cut-off LJ
with rcoulomb > rvdw for a while, but these were only available via
PME load balancing. Now we allow this setup to be chosen also
through mdp options.

Added optional support for portable hardware locality (hwloc)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added CMake support to detect and build |Gromacs| with hwloc, which
will improve |Gromacs| ability to recognize and take advantage of all
the available hardware. If hwloc is unavailable, |Gromacs| will fall back
on other detection routines.

Made normal-mode calculations work with shells and vsites
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Implemented shells and vsites in normal-mode analysis in mdrun and in
analysis of eigenvalues and frequencies. The normal-mode analysis
is done on real atoms only and the shells are minimized at each step
of the analysis.

:issue:`879`

Changed pull group count for coords stored in tpr file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added a parameter ngroup to the pull coord parameters. This is now
also stored in the tpr file. This makes the pull geometry forward
compatible, which is useful since it avoid bumping the .tpr version
with every new geometry, and we expect that users want to experiment
with new geometries.

Added pull coordinate geometry angle-axis
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The new geometry is described in the docs.
Some checks in readpull.cpp where reorganized since adding new
geometries made some old logic a bit convoluted.

Added pull coordinate geometry dihedral (angle)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
How to use the new geometry is explained in the docs.

Added pull coordinate geometry angle
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A new subsection was added to the docs explaining the new geometry.

Replaced ``pull-print-com1,2`` mdp option with ``pull-print-com``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Changes were made to the pull output order and naming.

Added pull potential flat-bottom-high
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added the new pull coordinate type flat-bottom-high, which is a flat
potential above the reference value and harmonic below.

Added ``gmx grompp`` check for pull group
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added a check for valid pull groups in a pull coordinate.
Using a pull group index that was out of range would cause invalid
memory access.

Added new swapping functionality to computational electrophysiology module
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Support was added for ion/water position swapping for multiple ion
types and polyatomic ions, including use of a user-defined number of
ionic species, and (small) polyatomic ions.

Also added two extra .mdp file parameters 'bulk-offset' that allow the
user to specify an offset of the swap layers from the compartment
midplanes. This is useful for setups where e.g. a transmembrane
protein extends far into at least one of the compartments. Without an
offset, ions would be swapped in the vicinity of the protein, which is
not wanted. Adding an extended water layer comes at the cost of
performance, which is not the case for the offset solution.

Documentation and testing was improved.

Fixed logic for DD missing-interactions check
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The code that was intended to double check that the domain decomposition
algorithm has not missed any interactions was inactive in several
cases, and has been fixed.

:issue:`1882`, :issue:`1793`

Permitted forces and velocities to be written to compressed TNG
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If there is no uncompressed coordinate output, write forces
and velocities to the TNG file with compressed coordinate
output. If there is uncompressed coordinate output to a
TNG file, forces and velocities will be written to it.

Use a greatest common divisor to set the frequency of some TNG
data output to ensure lambdas and box shape are written at least
as often as anything else.

:issue:`1863`

Added new notes to the user when coupling algorithms are unavailable
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
mdrun will now give the user an explanatory note when pressure and/or
temperature coupling is turned off.

Added mdrun check for finite energies
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added a check that the total potential energy is finite. This check is
nearly free and can catch issues with incorrectly set up systems
before users get a confusing constraint or PME error. Note that this
check is only performed at steps where energies are calculated, so it
will often not catch an exploding system.

Added ``gmx grompp`` check for unbound atoms
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
``gmx grompp`` now prints a note for atoms that are not connected by a
potential or constraint to any other atom in the same moleculetype,
since this often means the user made a mistake.

:issue:`1958`

Improved multi-simulation signalling
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Multi-simulations (including REMD) may have need to send messages
between the simulations. For example, REMD needs to write a
fully-consistent set of checkpoint files so that the restart works
correctly, but normal multi-simulations are fine with decoupled
progress and will simulate more efficiently if they can do
so. Similarly, ``gmx_mpi mdrun -maxh -multi`` needs to synchronize
only for REMD. The implementation has been upgraded so that such
coupling happens only when an algorithm chosen by the user requires
it.

:issue:`860`, :issue:`692`, :issue:`1857`, :issue:`1942`

Changed multi-simulation nsteps behaviour
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""-

It is unclear what the expected behaviour of a multi-simulation should
be if the user supplies any of the possible non-uniform distributions
of init_step and nsteps, sourced from any of .mdp, .cpt or command
line. Previously mdrun adjusted the total number of stesps to run so
that each run did the same number of steps, but now it reports on the
non-uniformity and proceed, assuming the user knows what they are
doing.

:issue:`1857`

Added working directory to things reported in .log file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
When running |Gromacs| via a batch script, it is useful to know which
working directory is being used for relative paths (file names) in the
command line. This is now written alongside other header information.

Prevented fragile use cases involving checkpoint restarts and/or appending
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

All output files named in the checkpoint file (ie. that were
used in the previous run) must be present before a checkpoint
restart will be permitted. Thus,
workflows where people used things like
``gmx mdrun -s production -cpi equilibration``
are no longer available to do a "continuous" restart. Instead, use
``gmx grompp -t equilibration -o production``.

:issue:`1777`

Removed warning after OpenMP core-count check
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
In many cases ``gmx_mpi mdrun`` issued a warning that compared the total
core count with something different returned from OpenMP. This problem
is caused by inappropriate management of thread affinity masks, but
the wording of the message did not help the user realise this, so has
been removed. ``gmx_mpi mdrun -pin on`` may help improve performance in
such cases.

Preparation for hardware detection might try to force offline cores to work
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Hardware detection might be foiled by kernels that take cores offline
when work is unavailable. We are not aware of any such platforms on which
|Gromacs| is likely to be used, but we will probably start to see them
soon. On such platforms, if the number of cores physically present
differs from the count that are online, we try to force them online
before deciding how |Gromacs| will use the online cores. For now, no x86 or
PowerPC platforms need such code, so it will never run on those platforms.
The good news is that we no longer have to risk making a confusing warning
about such possibilities.

Added new suggestion for users to try e.g. hyper-threading, if its disabled
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
|Gromacs| tends to perform best with several hardware threads available
per core (e.g. hyper-threading turned on, on x86), and now the log file
will note explicitly when such opportunities exist.


