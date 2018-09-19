Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

Implemented update groups
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Domain decomposition can now be based on so-called update groups. These
are groups of atoms with dependencies during the update, which can be
constraints and virtual sites. Update groups can typically be used when
only bonds involving hydrogens are constrained and are enabled
automatically when possible. This improves performance by eliminating
MPI and OpenMP communication for constraints and virtual sites.

PME on GPU when running free energy perturbations not involving charges
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
PME can now be run on a GPU when doing free energy perturbations
that do not involve perturbing charges.
