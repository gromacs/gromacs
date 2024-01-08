Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Improved Einstein viscosity calculation in gmx energy
"""""""""""""""""""""""""""""""""""""""""""""""""""""

Viscosity calcalution using the Einstein formula is convenient as this does
not require extremely frequent pressure tensor data. However, the implementation
of the calculation was inconvienent for long simulations and could take hours
to complete. Improved stepping through the data reduces the computational time
to minutes and provides much clearer output.

XVG output from ``gmx rdf`` now uses 6 decimal places
"""""""""""""""""""""""""""""""""""""""""""""""""""""

The output from ``gmx rdf`` now uses more decimal places in order to
avoid rounding issues. These issues led to perceived erroneous shifts in
the results.

:issue:`4647`

Handle CYX-CYX disulfide bonds in ``gmx pdb2gmx``
"""""""""""""""""""""""""""""""""""""""""""""""""

Naming CYS residues as CYX shows that they should form a disulfide bond.
``gmx pdb2gmx`` will now correctly interpret them as disulfide bond forming
residues.

:issue:`4929`

