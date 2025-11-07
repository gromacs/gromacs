New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


Collective variables module (Colvars) update
""""""""""""""""""""""""""""""""""""""""""""

The (`Colvars <https://colvars.github.io>`_) library for enhanced sampling simulations included
in |Gromacs| has been updated to version 2025-10-13.

This update brings many improvements, including:
- performance improvements for moving frame of reference calculations
- improved OPES implementation
- improved moving restraint logic, allowing for moving harmonic walls
- more flexible definitions of grids on CVs
A complete list of changes can be found `here <https://gitlab.com/gromacs/gromacs/-/merge_requests/5397>`_.

For documentation see :ref:`this section <colvars>`
of the |Gromacs| doc along with the `Colvars documentation page
<https://colvars.github.io/gromacs-2026/colvars-refman-gromacs.html>`_ dedicated to |Gromacs| 2026.
Messages about Colvars in the |Gromacs| discussion forum should be tagged
with the `colvars keyword <https://gromacs.bioexcel.eu/tag/colvars>`_ for
easier consultation.

Added LEaP-compatible dihedral reordering
"""""""""""""""""""""""""""""""""""""""""

AMBER LEaP reorders (improper) dihedrals when processing the topology, alphabetically
by atom type. In order for GROMACS to be able to match the topology produced by LEaP,
:ref:`gmx grompp` was extended LEaP-compatible dihedrals reordering. This functionality
is intended for usage with the ports of recent AMBER force fields (e.g. ff14SB and
ff19SB), but it is enabled on any force field that defines the
``_FF_AMBER_LEAP_ATOM_REORDERING`` macro.

This feature has been validated on all supported dipeptides and tripeptides as well as
select tetrapeptides and pentapeptides. Systematic validation on a broader set of
scientifically relevant biomolecules is still pending.

:issue:`4998`
