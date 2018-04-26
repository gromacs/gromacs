Topologies
==========

|Gromacs| must know on which atoms and combinations of atoms the various
contributions to the potential functions (see chapter :ref:`ff`) must act.
It must also know what parameters must be applied to the various
functions. All this is described in the *topology* file :ref:`top`, which
lists the *constant attributes* of each atom. There are many more atom
types than elements, but only atom types present in biological systems
are parameterized in the force field, plus some metals, ions and
silicon. The bonded and special interactions are determined by fixed
lists that are included in the topology file. Certain non-bonded
interactions must be excluded (first and second neighbors), as these are
already treated in bonded interactions. In addition, there are *dynamic
attributes* of atoms - their positions, velocities and forces. These do
not strictly belong to the molecular topology, and are stored in the
coordinate file :ref:`gro` (positions and velocities), or
trajectory file :ref:`trr` (positions, velocities, forces).

This chapter describes the setup of the topology file, the :ref:`top` file and
the database files: what the parameters stand for and how/where to
change them if needed. First, all file formats are explained. Section
:ref:`fffiles` describes the organization of the files in each force
field.

**Note:** if you construct your own topologies, we encourage you to
upload them to our topology archive at our `webpage`_! Just imagine how thankful
you’d have been if your topology had been available there before you
started. The same goes for new force fields or modified versions of the
standard force fields - contribute them to the force field archive!

.. _homepage: `webpage`_

.. toctree::
   :maxdepth: 2

   topologies/particle-type
   topologies/parameter-files
   topologies/molecule-definition
   topologies/constraint-algorithm-section
   topologies/pdb2gmx-input-files
   topologies/topology-file-formats
   topologies/force-field-organization


