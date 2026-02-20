Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Replaced usage of custom Bohr radius value in gmx spatial with the common value from units
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :ref:`gmx spatial` command used to have its own definition of Bohr radius. For consistency
with other parts of |Gromacs|, it now uses the definition of Bohr radius from the same source as
the rest of the code. Notably, the value of the constant changed from ``0.529177249`` (IUPAC 1999)
to ``0.529177210903`` (NIST 2018).

