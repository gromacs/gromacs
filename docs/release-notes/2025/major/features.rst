New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

A feature-limited version of the PLUMED interface is available
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A basic version of the `PLUMED <https://www.plumed.org/>`_ interface is now bundled by default in |Gromacs|.
With a non-Windows installation, it is possible to invoke PLUMED directly from the command line using the
``-plumed`` option of the ``gmx mdrun`` command, followed by the path to a PLUMED input file.
This can be done **without the need to apply a patch** as in previous |Gromacs| versions.
Importantly, this interface is not feature complete, see :ref:`the section in the manual <plumed>` for the details.

Support for amino-acid-specific energy correction maps (CMAPs)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Previously, energy correction map (CMAP) types could only be specified using
atom types and were applied to atoms in an amino-acid-agnostic way. CMAP types
can now be specified using both atom types and residue types, which enables
support for recent Amber force fields (ff19SB and later versions).

:issue:`4430`

Neural Network Potential support
""""""""""""""""""""""""""""""""
Basic support has been added to perform simulations with Neural Network Potentials (NNPs).
These models can be trained to reproduce forces and energies at *ab initio* levels of accuracy,
based on training data from electronic structure calculations with e.g. DFT or CCSD(T).
Importantly, |Gromacs| does *not* include any pretrained models, so users need to train their own
models or load pre-trained models from external sources.
As of now, the interface supports NNP models trained in `PyTorch <https://pytorch.org/>`_.
For details on usage and building |Gromacs| with LibTorch support, please see
the :ref:`NNPot section in the reference manual <nnpot>`.

Add Custom Improper Dihedrals in specbond.dat
"""""""""""""""""""""""""""""""""""""""""""""
This change allows users to specify an improper dihedral resulting for a special bond (i.e. a
thioester connection resulting in a SP2 group) in the specbond.dat file which is read during the
pdb2gmx preprocessing step.

:issue:`5113`
