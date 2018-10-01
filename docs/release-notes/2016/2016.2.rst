GROMACS 2016.2 Release Notes
----------------------------------------

This version was released on February 7, 2016. These release notes
document the changes that have taken place in GROMACS since version
2016.1 to fix known issues. It also incorporates all fixes made in
version 5.1.4 and several since.

Fixes where mdrun could behave incorrectly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add grompp check for equipartition violation risk for decoupled modes
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
When atoms involved in an angle with constrained bonds have very
different masses, there can be very weakly coupled dynamics modes.
Default mdp settings are often not sufficiently accurate to obtain
equipartitioning. This change adds a grompp check for this issue.

Part of :issue:`2071`

Disallow overwriting of dihedral type 9
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
It is no longer allowed to repeat blocks of parameter of multiple
lines for dihedrals of type 9. It is also no longer allowed to
override parameters or dihedrals of type 9. Both are too complex
to properly check. It is still allowed to repeat parameters blocks
consisting of a single line.
Repeating a block with the same parameters would lead to incorrect
dihedral potentials and forces.

:issue:`2077`

Fixed flat-bottom position restraints + DD + OpenMP
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A (re)allocation was missing, causing a crash.

:issue:`2095`

Fixed multi-domain reruns
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Old code cleanup led multi-domain rerun to crash because it failed to
consider logic separated over two places.

:issue:`2105`

Fixes for mdrun performance issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Corrected CUDA sm_60 performance
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The kernel launch now suits the SM size of the GP100 architecture.

Fixes for ``gmx`` tools
^^^^^^^^^^^^^^^^^^^^^^^

Fixed some FFT handling in cross-corrrelation calculations
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
An array of complex number was created as an array of pointers and
then passed to gmx_fft_1d. This does not work.

:issue:`2109`

Fixed gmx rmsf -q -oq
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This led to the PDB file containing B-factors using coordinates based
on those from the -s file, rather than -q file. gmx rmsf -oq was
otherwise fine.

Fixed crash in gmx order
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
gmx order used a cumbersome floating point method to compute
a histogram, leading to an index value that could be negative.

:issue:`2104`

Fixed minor trjconv bug
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
gmx trjconv -novel -f in.pdb -o out.pdb probably works better now.

Fixed time label print in gmx vanhove
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Handled issuing warnings correctly in xpm2ps and membed
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The code should not (over)write the output file before checking for
errors. For membed, it is useful to require the user to fix issues in
their input file before we unilaterally over-write it.

Corrected documentation about eigenvalue handling
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Some file format docs were out of step with the implementation in
eigio.cpp.

The behaviour of gmx anaeig -eig -eig2 was not properly documented.

Made editconf B-factor attachment more useful in practice
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
B-factor values will be added to residues unless an index is larger
than the number of residues or an option is specified. Protein residue
indices can start from any number and, in case they start from a large
number, there is no way to add B-factor values to residues.

This patch changes it to add B-factor values to residues unless the
number of B-factor values is larger than the number of residues.

Fixed possible memory error with long selections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If a selection was more than 1000 characters long and there was a
whitespace exactly at the 1000 point, a buffer overflow could occur.
Replaced the buffer with std::string, simplifying the code
significantly.

:issue:`2086`

Fixed use of position variables with plus/merge
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If a selection contained a position variable (e.g., 'com of ...') that
was used more than once, and at least one of those uses was with
plus/merge, there were out-of-bounds memory writes.  This was caused by
the internal position structure not getting fully initialized.
Incomplete initialization happens in all contexts with such variables,
but only plus/merge (and possibly permute) actually use the values that
remained uninitialized, which caused them to incorrectly compute the
amount of memory required to store the result.

:issue:`2086`

Improved documentation
^^^^^^^^^^^^^^^^^^^^^^

Made several minor improvements to documentation and messages to users
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
In particular, for selections:

- Explained resindex and resnr keywords in selection help.
- Explained how selection-enabled tools treat -s and -f input files.

:issue:`2083`

Clarified use of tau-p and pcoupltype
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
grompp used to permit the erroneous "tau-p = 5 5". This does not
reflect that only one time constant is permitted for pressure coupling
(unlike group-based temperature coupling). The recent fix for
:issue:`1893` leads to the user receiving a grompp warning, so this
improves the docs to make clear that pressure coupling is different.

:issue:`1893`

Portability enhancements
^^^^^^^^^^^^^^^^^^^^^^^^

Fixed x86 conditional on IBM s390x
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The CpuInfoTest.SupportLevel test fails on IBM s390x because wrong
condition was used.

Fixes: https://bugzilla.redhat.com/show_bug.cgi?id=1390149

:issue:`2072`

Build system enhancements
^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed compilation with CMAKE_CXX_FLAGS="-Wall -Werror"
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`2073`

Stopped trying to use objdump --reloc in the build system on Mac
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Recent Xcode objdump does not support --reloc.

The warning that is based on the output of running objdump was only
implemented to work on Linux-like things, so should not spam the cmake
output on other platforms.

Improved the support for plugin loading in the build system
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The mdrun-only and prefer-static-libs builds set the default for
BUILD_SHARED_LIBS to off, which silently disabled plugin support
for things like VMD-based I/O handling.

Converted GMX_LOAD_PLUGINS to tri-state ON/OFF/AUTO so that if the
preconditions for support are not met we can have suitable behaviour
in each case.

:issue:`2082`

Turn off hwloc support when static lib found
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Hwloc dependencies are not resolved at CMake time when static
libwloc.a is detected and in most of these cases link-time
errors will prevent building GROMACS. As it is hard for a user to know
how to solve such cryptic errors and hwloc is not a required dependency,
we turn off hwloc support when a static lib is detected. The user can
override this on the cmake command line.

:issue:`1919`

Fixed build with GMX_EXTERNAL_TNG=ON
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

House-keeping that reduces users' problems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mdrun prints invalid performance data less often
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
If mdrun finished before a scheduled reset of the timing information
(e.g. from mdrun -resetstep or mdrun -resethway), then misleading
timing information should not be reported.

Related, the default reset step for gmx tune_pme was increased to 1500.

:issue:`2041`

Added a runtime check for number of threads in bonded code
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Replaced a debug assertion on the number of OpenMP threads not being
larger than GMX_OPENMP_MAX_THREADS by fatal error.
But since the listed-forces reduction is actually not required with
listed forces, these are now conditional and mdrun can run with more
than GMX_OPENMP_MAX_THREADS threads.

:issue:`2085`

Fixed integer narrowing in TNG reading for long trajectories
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Reading of TNG trajectories with sufficiently large numbers of frames
could truncate integers used for frame numbers. Fixed to use 64-bit
integers as originally intended.

Fixed logic of TRR reading
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
When reading a trr file, reaching the end of the file was
indistinguishable from a reading error or a magic-number error. This
is now fixed, restoring the intended behaviour in each case.

:issue:`1926`
