Removed mdrun features
^^^^^^^^^^^^^^^^^^^^^^

Removed SD2 integrator
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This integrator has known problems, and is in all ways inferior to
sd. It has no tests, and was deprecated in |Gromacs| 5.0. There are no
plans to replace it.

:issue:`1137`

Removed the twin-range scheme
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Only the (deprecated) group scheme supports this, and the Verlet scheme will not
support it in the foreseeable future.  There is now the explicit
requirement that rlist >= max(rcoulomb,rvdw).

Removed support for twin-range with VV integrators
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Group-scheme twin-ranged non-bonded interactions never worked with
velocity-Verlet integrators and constraints. There are no plans to
make that combination work.

:issue:`1137`, :issue:`1793`

Removed Reaction-Field-nec
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The reaction-field no-exclusion correction option was only introduced for
backward compatibility and a performance advantage for systems
with only rigid molecules (e.g. water). For all other systems
the forces are incorrect. The Verlet scheme does not support this
option and even if it would, it wouldn't even improve performance.

Removed AdResS module
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This feature requires the (deprecated) group scheme, and there are no
plans to port it to the Verlet scheme.

:issue:`1852`

Removed mdrun -compact
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
It is too complicated to support multiple ways of analysing per-step
data.

Removed lambda printing from mdrun log file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`1773`

Removed GMX_NOCHARGEGROUPS
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This undocumented feature was only useful with the (deprecated) group
scheme.
