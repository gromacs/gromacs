New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

Update gmx cluster to write correct PDB files and index files with cluster frames
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:ref:`PDB <pdb>` files from gmx cluster were missing the CRYST header for box information, making
it more difficult than needed to use them with our |Gromacs| tools. Also, the :ref:`index <ndx>`
files needed for :ref:`gmx trjconv` to split up trajectories into frames corresponding
to the clusters were not written. This adds support for writing this :ref:`index <ndx>` file
as well as proper :ref:`PDB <pdb>` files.
