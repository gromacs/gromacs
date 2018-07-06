Miscellaneous
^^^^^^^^^^^^^

pdb2gmx writes total charge differently
---------------------------------------

In the ``[atoms]`` field of the topology files it writes, pdb2gmx
writes the total charge for each residue. This should generally be an
integer, so can help with troubleshooting issues in system or force
field preparation. That printing is now done only once per residue,
rather than for every atom.
