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
``gmx x2top`` can write MINFF topologies and has new options
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

:ref:`gmx x2top` has been extended to make topologies for mineral systems with the
`MINFF <https://github.com/mholmboe/minff>`_ force field (Holmboe, Lavauzelle and Swai,
arXiv preprint, 2026, doi:10.48550/arXiv.2610.18220). Some of the new options are
also useful with other force fields, and none of them changes the output unless it is given.

* ``-n2t`` reads the atom name to atom type table from a given file, instead of
  from the ``atomname2type.n2t`` files of the force field.
* ``-ff min`` (or ``-ff minff``) writes the compact MINFF format: bonds to
  hydrogen only, and angles with the measured angle and force constant, annotated with the
  atom types. It sets the defaults of ``-kb``, ``-kt``, ``-ktH``, ``-nexcl`` and ``-name``
  for MINFF. The MINFF force field directory ``min.ff`` is not distributed with |Gromacs|.
* ``-ktH`` sets a separate angle force constant for angles that involve a hydrogen
  atom. Hydrogens are recognized from the mass given in the n2t table, otherwise from
  the element column of a PDB file, and otherwise from the atom name.
* ``-dH``, ``-dM`` and ``-aH`` override the printed bond length of bonds to hydrogen, the
  bond length of all bonds, and the angle of angles involving hydrogen, in the MINFF format.
* ``-smear`` (with ``-ff min``) computes the charges of oxygens and structural fluorines
  from the charges of their neighbors, instead of reading them from the n2t table. This
  smears out the charge of isomorphic substitutions, see Tournassat et al.,
  Clays and Clay Minerals 64, 374-388 (2016), doi:10.1346/CCMN.2016.0640403 and
  Lammers et al., J. Colloid Interface Sci. 2017, doi:10.1016/j.jcis.2016.11.084.
  The total charge is rounded to the nearest integer.

Bonds, angles and dihedrals are now found using the periodic boundary condition type
of the input structure, instead of always assuming full periodicity for angles and
dihedrals. Output for fully periodic input is unchanged.
