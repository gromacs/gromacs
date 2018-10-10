Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pdb2gmx writes total charge differently
---------------------------------------

pdb2gmx notes the total charge for each residue in the ``[atoms]``
field of the topology file it produces. The fact that this should
generally be an integer can be used for troubleshooting issues in
system or force field preparation. This printing is now done only once
per residue, rather than for every atom.

nmeig does thermochemistry
---------------------------------------

The nmeig tool that analyzes the Hessian matrix from a normal mode
analysis now generates thermochemical properties like standard
entropy, heat capacity at constant volume, internal thermal energy
and zero-point energy. The analysis is based on the harmonic
approximation that is the same as what is used in quantum chemistry.

Implement writing of LaTeX methods in report-methods
----------------------------------------------------

Added a new tool report-methods to write
a summary of methods used to set up a simulation to an output file
or to standard output.
