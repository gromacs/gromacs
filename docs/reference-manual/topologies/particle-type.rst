Particle type
-------------

In |Gromacs|, there are three types of
particles
, see :numref:`Table %s <tab-ptype>`. Only regular atoms and virtual
interaction sites are used in |Gromacs|; shells are necessary for
polarizable models like the Shell-Water models \ :ref:`45 <refMaaren2001a>`.

.. _tab-ptype:

.. table:: Particle types in |Gromacs|

           +--------------+----------+
           | Particle     | Symbol   |
           +==============+==========+
           | atom         | A        |
           +--------------+----------+
           | shell        | S        |
           +--------------+----------+
           | virtual side | V (or D) |
           +--------------+----------+


.. _atomtype:

Atom types
~~~~~~~~~~

Each force field defines a set of atom
types,
which have a characteristic name or number, and mass (in a.m.u.). These
listings are found in the ``atomtypes.atp`` file (:ref:`atp` =
**a**\ tom **t**\ ype **p**\ arameter file). Therefore, it is in this
file that you can begin to change and/or add an atom type. A sample from
the ``gromos43a1.ff`` force field is listed below.

::

     |  O  15.99940 ;     carbonyl oxygen (C=O)
     | OM  15.99940 ;     carboxyl oxygen (CO-)
     | OA  15.99940 ;     hydroxyl, sugar or ester oxygen
     | OW  15.99940 ;     water oxygen
     |  N  14.00670 ;     peptide nitrogen (N or NH)
     | NT  14.00670 ;     terminal nitrogen (NH2)
     | NL  14.00670 ;     terminal nitrogen (NH3)
     | NR  14.00670 ;     aromatic nitrogen
     | NZ  14.00670 ;     Arg NH (NH2)
     | NE  14.00670 ;     Arg NE (NH)
     |  C  12.01100 ;     bare carbon
     |CH1  13.01900 ;     aliphatic or sugar CH-group
     |CH2  14.02700 ;     aliphatic or sugar CH2-group
     |CH3  15.03500 ;     aliphatic CH3-group

**Note:** |Gromacs| makes use of the atom types as a name, *not* as a
number (as *e.g.* in GROMOS).

.. _vsitetop:

Virtual sites
~~~~~~~~~~~~~

Some force fields use virtual interaction sites (interaction sites that
are constructed from other particle positions) on which certain
interactions are located (*e.g.* on benzene rings, to reproduce the
correct quadrupole). This is described in sec. :ref:`virtualsites`.

To make virtual sites in your system, you should include a section
``[ virtual_sites? ]`` (for backward compatibility the old
name ``[ dummies? ]`` can also be used) in your topology
file, where the ``?`` stands for the number constructing
particles for the virtual site. This will be ``2`` for
type 2, ``3`` for types 3, 3fd, 3fad and 3out and
``4`` for type 4fdn. The last of these replace an older
4fd type (with the ‘type’ value 1) that could occasionally be unstable;
while it is still supported internally in the code, the old 4fd type
should not be used in new input files. The different types are explained
in sec. :ref:`virtualsites`.

Parameters for type 2 should look like this:

::

    [ virtual_sites2 ]
    ; Site  from        funct  a
    5       1     2     1      0.7439756

for type 3 like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          b
    5       1     2     3      1       0.7439756  0.128012

for type 3fd like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          d
    5       1     2     3      2       0.5        -0.105

for type 3fad like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   theta      d
    5       1     2     3      3       120        0.5

for type 3out like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          b          c
    5       1     2     3      4       -0.4       -0.4       6.9281

for type 4fdn like this:

::

    [ virtual_sites4 ]
    ; Site  from                      funct   a          b          c
    5       1     2     3     4       2       1.0        0.9       0.105

This will result in the construction of a virtual site, number 5 (first
column ``Site``), based on the positions of the atoms
whose indices are 1 and 2 or 1, 2 and 3 or 1, 2, 3 and 4 (next two,
three or four columns ``from``) following the rules
determined by the function number (next column ``funct``)
with the parameters specified (last one, two or three columns
``a b . .``). Obviously, the atom numbers (including
virtual site number) depend on the molecule. It may be instructive to
study the topologies for TIP4P or TIP5P water models that are included
with the |Gromacs| distribution.

**Note** that if any constant bonded interactions are defined between
virtual sites and/or normal atoms, they will be removed by
:ref:`grompp <gmx grompp>` (unless the option ``-normvsbds`` is used). This
removal of bonded interactions is done after generating exclusions, as
the generation of exclusions is based on “chemically” bonded
interactions.

Virtual sites can be constructed in a more generic way using basic
geometric parameters. The directive that can be used is ``[ virtual_sitesn ]``. Required
parameters are listed in :numref:`Table %s <tab-topfile2>`. An example entry for
defining a virtual site at the center of geometry of a given set of
atoms might be:

::

    [ virtual_sitesn ]
    ; Site   funct    from
    5        1        1     2     3     4
