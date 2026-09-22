New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Multiple molecule types with identical SETTLE parameters are now supported
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This can be useful when one, for instance, wants to apply position restraints
to a subset of water molecules.

Periodic SCC-DFTB, GFN1-xTB and GFN2-xTB are available for CP2K QM/MM
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The CP2K QM/MM interface can now generate periodic SCC-DFTB, GFN1-xTB and
GFN2-xTB inputs using short-range regularized point-charge coupling with SPME
electrostatics. For GFN1-xTB and GFN2-xTB CP2K should be compiled with support 
for the tblite interface. 
The user can select Gaussian Expansion of Electrostatic Potential (GEEP) 
for MM charges in tight-binding methods through the
``qmmm-cp2k-dftb-electrostatic-coupling = gauss`` MDP option, 
which requires a CP2K 2027.1 or higher, whereas the default is 
short-range regularized point-charge coupling 
``qmmm-cp2k-dftb-electrostatic-coupling = point-charge``.
Regular DFT methods are still using GEEP (Gauss) for electrostatics.

CMAP interactions now support free energy perturbation
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Energy correction map (CMAP) torsion interactions can now be
perturbed between an A state and a B state using the
``bonded-lambdas`` free-energy lambda component. The A and B
state CMAP grids are linearly interpolated at each lambda
value, and the corresponding :math:`\partial H/\partial\lambda`
contribution is accumulated for BAR/TI analysis.

To use CMAP FEP, specify a B-state CMAP type in the ``[ cmap ]``
section of the topology. The type can be given as a 1-based index
or as an explicit name token (e.g. ``GLY``, ``PZQ``) that matches
the optional name field on the corresponding ``[ cmaptypes ]``
header line. When the A and B state CMAP types differ, the
interaction is treated as perturbed; when they are identical (or
no B-state is specified) the standard non-perturbed code path is
used.
