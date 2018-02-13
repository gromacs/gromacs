Miscellaneous
^^^^^^^^^^^^^

Updated note in manual on stochastic dynamics integrator
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The comment in the SD section about Berendsen was outdated.
Added a few sentences on equilibration/damping of modes.

Added grompp note for Parrinello-Rahman + position restraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This combination can be unstable and is often not desirable, so
grompp now issues a note to suggest alternatives to the user.

Refs :issue:`2330`

Clarified the description of Fmax during energy minimization
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Improved vsite parallel checking
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The vsite struct now stores internally whether it has been configured
with domain decomposition. This allows for internal checks on valid
commrec, which have now been added, and would have prevented :issue:`2257`.

Added partial support for writing masses and partial charges with TNG files
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:issue:`2188`

Updated TNG to version 1.8.1
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Added data block for atom masses.

Fixes :issue:`2187` and :issue:`2250` and other bugs and warnings.

Added load balance fraction to DLB print
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
DLB can often be based on a small fraction of the total step time,
especially with GPUs. Now this is printed to md.log and stderr.

Added reference for dihedral function in OPLS.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The OPLS four-term dihedral function was not described in the
reference listed earlier, so this was updated. Also updated
the reference to the three term dihedral to an older paper.

Updated many aspects of the documentation
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Imported and updated more material from the wiki. Incorporated
suggestions arising from many Redmine issues. Updated user guide,
developer guide, install guide, and reference manual.

Updated mdrun signal help text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Updated mdrun help text on signal handling for old and recent changes
to the behavior.

Fixes :issue:`2324`

Changed to handle erroneous command line args better
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Some gmx modules need to be able to accept non-option arguments, and
some should not. Introduced enough functionality to support such
behaviour, while giving useful error messages in cases where the
command line is merely missing hyphens (which can happen e.g. when
people copy-paste from inconveniently built PDF files for tutorials).
Increased test coverage of relevant cases.

Removed some useless command-line argument strings from test cases
that never needed them.

Also tested some behaviours of handling string options, and renamed
some test input strings to reflect the intent.

:issue:`2153`

Changed to no longer allow multiple energy groups for GPU runs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Exit with a fatal error instead of only warning, since the latter
leads to writing data for energy groups that is incorrect to the
energy file.

:issue:`1822`

Removed duplications in GMXLIB search paths
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Remove entries that are duplicated, or identical to the default search
path, to avoid e.g.  listing identical force fields multiple times.

:issue:`1928`

Changed to no longer write reference pull group 0 to log
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This is an internal group used for absolute references, which cannot
be set by users, so printing it just leads to confusion.

:issue:`2143`
