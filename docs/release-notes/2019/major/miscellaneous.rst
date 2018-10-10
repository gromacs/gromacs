Miscellaneous
^^^^^^^^^^^^^

grompp discourages use of constraints=all-bonds
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Common force fields, including AMBER, CHARMM and OPLS/aa, are parametrized
with bonds involving hydrogen constrained. Constraining all bonds should
be avoided, for correctness. grompp now issues a note when
constraints=all-bonds is used with these force fields when time steps
are smaller than 2.6 fs and hydrogens are not replaced by virtual sites.
Using constraints=h-bonds will also improve performance.

Documentation changed to reStructuredText
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The complete documentation has been moved to the reStructuredText 
markup format to allow building it together for display as html or as pdf.
