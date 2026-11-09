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

``-w`` now displays more file types, and falls back to other viewers
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The ``-w`` option can now display ``.cif`` and ``.gro`` files, neither of which
had a viewer before. Every file type it handles now lists several viewers rather
than one: the viewer used previously is still preferred, and an alternative is
chosen only when that one is not installed, so ``-w`` now works on machines where
it used to fail. Structures are shown with PyMOL, ChimeraX, Jmol or VMD, plots
with Grace, and images and PostScript with the common desktop viewers. The
``GMX_VIEW_`` environment variables continue to override the choice.

``gmx trjconv`` now supports converting trajectories to H5MD
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

H5MD is now supported as both an input and output format for trajectory
conversion.
