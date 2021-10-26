Bugs fixed
^^^^^^^^^^

Fixed slight inaccuracies when using virtual sites with pressure coupling
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Virtual sites were reconstructed after the system was propagated, but before
scaling due to pressure coupling. For virtual site types which are not a linear
combination of other atoms, this is not completely correct. Since the scaling
due to pressure coupling is very small in healthy simulations, the resulting
inaccuracies are expected to have been extremely minor, and in most cases
undetectable.

:issue:`3866`

Correct dVremain/dl when nstdhdl > nstcalcenergy
""""""""""""""""""""""""""""""""""""""""""""""""

When nstcalcenergy was not a multiple of nstdhdl, incorrect dVremain/dl
terms were written in the energy file. Note that all dH/dl output in
both dhdl.xvg and the energy file, which is used by e.g. gmx bar, was correct.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Use correct c0 parameter in Me2PO4 in OPLSAA
""""""""""""""""""""""""""""""""""""""""""""

OPLSAA torsions must sum to 0, but the paramters for Me2PO4 did not do so. Changed the c0
parameter to the correct value.

:issue:`4075`

