Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Extend supported use-cases for GPU version of update and constraints
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

GPU version of update and constraints can now be used for FEP, except mass and constraints
free-energy perturbation.
       
Reduce time spent in grompp with large numbers of distance restraints
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The time `gmx grompp` spent processing distance restraint has been
changed from quadratic in the number of restraints to linear.
       
:issue:`3457`

Support for offloading PME to GPU when doing Coulomb FEP
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

PME calculations can be offloaded to GPU when doing Coulomb free-energy perturbations.

CPU SIMD accelerated implementation of harmonic bonds
"""""""""""""""""""""""""""""""""""""""""""""""""""""

SIMD acceleration for bonds slightly improves performance for systems
with H-bonds only constrained or no constraints. This gives a significant
improvement with multiple time stepping.
