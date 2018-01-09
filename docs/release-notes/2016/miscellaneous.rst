Miscellaneous
^^^^^^^^^^^^^

Various improvements to documentation and tests
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In particular, the definition of pressure in the reference manual
should be in bar, and a spurious r_ij in the force for the Morse
potential was removed. Added documentation and literature references
for membrane embedding. Improved template analysis program
documentation. gmock was patched to work with gcc 6.

:issue:`1932`

Improved make_ndx help text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Clarified the use of boolean operators. The old help text could
incorrectly hint that AND, OR, and NOT would work as keywords.
Added a reference to ``gmx select`` that in most cases can serve as a
replacement.

:issue:`1976`

Addded checks on number of items read in mdp statements
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added checks for the number of items read in all
sscanf() statements processing data from the mdp
file.

:issue:`1945`.

Work around glibc 2.23 with CUDA
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
glibc 2.23 changed the behaviour of string.h in a way that broke all
versions of CUDA with all gcc compiler versions. The GROMACS build
system detects this glibc, and works around it by adding the
_FORCE_INLINE preprocessor define to CUDA compilation.

:issue:`1982`

Split NBNXN CUDA kernels into four compilation units
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The CUDA nonbonded kernels are now built in four different compilation units
when this is possible; ie. devices with compute capability >= 3.0. This
can dramatically reduce compilation time.

Forcing the use of a single compilation unit can be done using the
GMX_CUDA_NB_SINGLE_COMPILATION_UNIT cmake option.

:issue:`1444`

Added stream flushes when not writing newline character
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Some of our routines use the carriage return without a newline
to keep writing the status e.g. on stderr.
For some operating systems this seems to lead to the output
being cached in the buffers, so this change adds an explicit
fflush() for these print stamements.

Fixed :issue:`1772`

Supported cmap with QMMM
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Formerly, QMMM only supported bonded interactions using up to 4 atoms.
Now any number is supported and some hard-coded assumptions have been
removed.

Upgraded support for lmfit library
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Now based on lmfit 6.1. The CMake option GMX_EXTERNAL_LMFIT permits
linking an external lmfit package, rather than the one bundled in
GROMACS.

:issue:`1957`

libxml2 is no longer a dependency
""""""""""""""""""""""""""""""""""""""""""""""""""""""
GROMACS used to use libxml2 for running its test code. This has been
replaced by a bundled version of tinyxml2 (or optionally, a system
version of that library).

Disable automated FFTW3 builds on Windows
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The FFTW distribution does not include configurations to
build it automatically on windows, in particular not through
the ``./configure; make; make install`` triad.

:issue:`1961`

Remove warnings on checkpoint mismatch
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
mdrun now only warns for mismatch in minor version, build or
number of ranks used when reproducibility is requested.
Also added a separate message for not matching precision.

:issue:`1992`

Report the filename and the line number on failure
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Extend the call to gmx_fatal in fget_lines() to report the filename and
the line number where the read failed.

Handled constraint errors with EM
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
All energy minimizers could fail with random errors when constraining
produced NaN coordinates.
Steepest descents now rejects steps with a constraint error.
All other minimizer produce a fatal error with the suggestion to use
steepest descents first.

:issue:`1955`

Disable static libcudart on OS X
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Recent versions of CMake enable a static version of
libcudart by default, but this breaks builds at least
on the most recent version (10.11) of OS X, so we
disable it on this platform.

Fixed rare issue linking with clock_gettime
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Misuse of preprocessing commands might have led to inappropriate
use of clock_gettime().

:issue:`1980`

Disabled NVIDIA JIT cache with OpenCL
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The NVIDIA JIT caching is known to be broken with OpenCL compilation in
the case when the kernel source changes but the path does not change
(e.g. kernels get overwritten by a new installation). Therefore we disable
the JIT caching when running on NVIDIA GPUs. AMD GPUs are unaffected.


:issue:`1938`

