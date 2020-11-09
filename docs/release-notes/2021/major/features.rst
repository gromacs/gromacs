New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Virtual site with single constructing atom
""""""""""""""""""""""""""""""""""""""""""

Added a virtual site that is constructed on top if its single constructing
atom. This can be useful for free-energy calculations.

Density-guided simulations can apply matrix multiplication and shift vector to structures
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The new mdp option "density-guided-simulation-shift-vector" defines a
shift vector that shifts the density-guided simulation group before the 
density forces are evaluated. With a known shift vector that aligns structure
and input density, this feature enables structure refinement to non-aligned
densities without the need to manipulate the input density data or structure.
The mdp option "density-guided-simulation-transformation-matrix" allows to 
define a matrix with which to multiply the structure coordinates, before the shift
vector is applied. This allows arbitrary rotation, skewing and scaling of input
structures with respect to the input densities.
A typical use case are membrane-embedded proteins which cannot easily be
shifted and rotated within membranes.

Lower energy drift due to SETTLE
""""""""""""""""""""""""""""""""

|Gromacs| already applied an improvement to the center of mass calculation in
SETTLE to reduce energy drift in single precision. Now the center of mass
calculation is completely avoided, which significantly reduces the energy
drift when large coordinate values are present. This allows for accurate
simulations of systems with SETTLE up to 1000 nm in size (but note that
constraining with LINCS and SHAKE still introduces significant drift,
which limits the system size to 100 to 200 nm).

mdrun now reports energy drift
""""""""""""""""""""""""""""""

With conservative integrators, mdrun now reports the drift of the conserved
energy quantity in the log file.

FEP using AWH
"""""""""""""

It is now possible to control the lambda state of a free energy perturbation
simulation using the Accelerated Weight Histogram method. This can be used
as one of multiple AWH dimensions, where the other(s) are coupled to pull
coordinates.

Support for cyclic molecules in pdb2gmx
"""""""""""""""""""""""""""""""""""""""

It is now possible to process cyclic molecules in pdb2gmx and generate |Gromacs|
topology files for them.

Stochastic cell rescaling barostat
""""""""""""""""""""""""""""""""""

Implementation of the stochastic cell rescaling barostat. This is a first-order,
stochastic barostat, that can be used both for equilibration and production.
