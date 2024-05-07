Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Verlet buffer set correctly for inhomogeneous systems
"""""""""""""""""""""""""""""""""""""""""""""""""""""

The Verlet buffer estimation now uses an effective density for
the system computing from the initial coordinates. This avoids
underestimation of the buffer for (very) inhomogeneous systems.

:issue:`4509`

Fix segmentation fault for large atom and thread count
""""""""""""""""""""""""""""""""""""""""""""""""""""""

When the number of atoms times the number of OpenMP threads was larger
than 2147483647, negative atom number could cause illegal memory access.

:issue:`4628`

Density-guided simulation normalization
"""""""""""""""""""""""""""""""""""""""

With the ``.mdp`` option ``density-guided-simulation-normalize-densities = yes``
, the reference density and the simulated density values were previously divided
by the sum of their values.

This lead to surprising behavior for reference densities with lots of negative
voxel values: the density started to "repel" the protein structure
instead of attracting it, if the total sum of voxel values was smaller
than zero. The negative normalization constant lead to a sign change in voxel
values.

To avoid this behavior, the reference density is now normalized so that the
sum of *positive* values is unity, ensuring that the normalization constant is
always positive.

Apart from avoiding the unexpected behavior, we expect that this also leads 
to smaller absolute differences between reference density and simulated density,
with some small benefits for numerical stability.

This change affects all simulations where voxel values are negative
(usually this excludes synthetic data) and that are run with
``density-guided-simulation-normalize-densities = yes``, but only has a larger
effect for: first, similarity  measure ``inner-product`` as an effective
force-constant scaling and, second, for all similarity measures where the sum
of all voxel values was negative.   

gmxapi Python package avoids unnecessary MPI initialization
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Delayed initialization of MPI (due to automatic behavior of :py:mod:`mpi4py`)
avoids MPI initialization that previously occurred just by importing :py:mod:`gmxapi`.
The previous behavior has been seen to cause strange interactions with
resource management libraries like ``libfabric`` at unexpected times
(such as during package installation) with :py:mod:`gmxapi` version 0.3.

:issue:`4693`

Fail-safe check for perturbed exclusions beyond rlist
"""""""""""""""""""""""""""""""""""""""""""""""""""""

With free-energy calculations, excluded non-bonded interactions involving
at least one perturbed atom should not be beyond rlist when using PME. The
check for this could have false negatives. Now the check is fail safe and
will always trigger a fatal error when perturbed excluded pairs are beyond rlist.

:issue:`3403`
:issue:`4321`
:issue:`4461`
