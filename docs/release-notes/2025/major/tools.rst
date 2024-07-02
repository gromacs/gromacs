Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

``gmx grompp`` now checks dihedral coefficients sum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sum of dihedral parameters of type 3 (i.e Ryckaert-Bellemans or Fourier
dihedrals) is now checked during preprocessing. In free energy simulations,
this sum must be equal in both states as it affects final results 
through dH/dl. Additionally, this sum should be zero when comparing potential
energy values with other force field ports and simulation codes, but a non-zero
sum does not otherwise affect the simulation  (a simple note is emitted).

If ``gmx grompp`` rejects an interaction that was previously accepted,
then change the first coefficient to make the total zero. This leading
coefficient has no effect on the derivative of the energy (ie. the forces),
and thus no effect on the dynamics.

No parameters in default force fields in |Gromacs| were affected, so
none have changed.

:issue:`4253`

Added support for DSSP v.2 in ``gmx dssp``
""""""""""""""""""""""""""""""""""""""""""

There is now the ability to choose between two different output modes: with polyproline helices
search enabled (option "-polypro", default and corresponds to the output of DSSP v.4) and
with polyproline helices search disabled (option "-nopolypro", corresponds to the output of DSSP v.2).
