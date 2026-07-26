Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

``gmx pdb2gmx`` now resolves generated hydrogen names from the residue topology
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When ``gmx pdb2gmx`` adds hydrogens from the hydrogen database, it now looks up
their names in the residue topology entry of the selected force field instead of
relying on hard-coded numbering rules. This allows force fields that use
different hydrogen naming schemes, such as original AMBER/IUPAC names for
methylene hydrogens, to work without naming-specific code paths or command-line
options.

