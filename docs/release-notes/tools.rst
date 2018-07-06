Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pdb2gmx writes total charge differently
---------------------------------------

pdb2gmx notes the total charge for each residue in the ``[atoms]``
field of the topology file it produces. The fact that this should
generally be an integer can be used for troubleshooting issues in
system or force field preparation. This printing is now done only once
per residue, rather than for every atom.

