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
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""^

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

Updated developer guide
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Imported and updated more material from the wiki. Included coverage of
some recent discussion points on C++11 and preprocessor use.

Updated mdrun signal help text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Updated mdrun help text on signal handling for old and recent changes
to the behavior.

Fixes :issue:`2324`

