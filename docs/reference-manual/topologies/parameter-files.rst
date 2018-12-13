Parameter files
---------------

Atoms
~~~~~

The *static* properties (see  :numref:`Table %s <tab-statprop>`) assigned to the atom
types are assigned based on data in several places. The mass is listed
in ``atomtypes.atp`` (see :ref:`atomtype`), whereas the charge is listed
in :ref:`rtp` (:ref:`rtp` = **r**\ esidue **t**\ opology **p**\ arameter file,
see :ref:`rtp`). This implies that the charges are only defined in the
building blocks of amino acids, nucleic acids or otherwise, as defined
by the user. When generating a :ref:`topology <top>` using the :ref:`pdb2gmx <gmx pdb2gmx>`
program, the information from these files is combined.

.. _tab-statprop:

.. table:: Static atom type properties in |Gromacs|

           +----------+------------------+----------+
           | Property | Symbol           | Unit     |
           +==========+==================+==========+
           | Type     | -                | -        |
           +----------+------------------+----------+
           | Mass     | m                | a.m.u.   |
           +----------+------------------+----------+
           | Charge   | q                | electron |
           +----------+------------------+----------+
           | epsilon  | :math:`\epsilon` | kJ/mol   |
           +----------+------------------+----------+
           | sigma    | :math:`\sigma`   | nm       |
           +----------+------------------+----------+


.. _nbpar:

Non-bonded parameters
~~~~~~~~~~~~~~~~~~~~~

The non-bonded parameters consist of the van der Waals parameters V (``c6``
or :math:`\sigma`, depending on the combination rule) and W (``c12`` or
:math:`\epsilon`), as listed in the file ``ffnonbonded.itp``, where ``ptype`` is
the particle type (see :numref:`Table %s <tab-ptype>`). As with the bonded
parameters, entries in ``[ *type ]`` directives are applied to their counterparts in
the topology file. Missing parameters generate warnings, except as noted
below in section :ref:`pairinteractions`.

::

    [ atomtypes ]
    ;name   at.num      mass      charge   ptype         V(c6)        W(c12)
        O        8  15.99940       0.000       A   0.22617E-02   0.74158E-06
       OM        8  15.99940       0.000       A   0.22617E-02   0.74158E-06
       .....

    [ nonbond_params ]
      ; i    j func       V(c6)        W(c12)
        O    O    1 0.22617E-02   0.74158E-06
        O   OA    1 0.22617E-02   0.13807E-05
        .....

**Note** that most of the included force fields also include the ``at.num.``
column, but this same information is implied in the OPLS-AA ``bond_type``
column. The interpretation of the parameters V and W depends on the
combination rule that was chosen in the ``[ defaults ]`` section of the topology file
(see :ref:`topfile`):

.. math:: \begin{aligned}
          \mbox{for combination rule 1}: & &
          \begin{array}{llllll}
            \mbox{V}_{ii} & = & C^{(6)}_{i}  & = & 4\,\epsilon_i\sigma_i^{6} &
            \mbox{[ kJ mol$^{-1}$ nm$^{6}$ ]}\\
            \mbox{W}_{ii} & = & C^{(12)}_{i} & = & 4\,\epsilon_i\sigma_i^{12} &
            \mbox{[ kJ mol$^{-1}$ nm$^{12}$ ]}\\
          \end{array}
          \\
          \mbox{for combination rules 2 and 3}: & &
          \begin{array}{llll}
            \mbox{V}_{ii} & = & \sigma_i   & \mbox{[ nm ]} \\
            \mbox{W}_{ii} & = & \epsilon_i & \mbox{[ kJ mol$^{-1}$ ]}
          \end{array}\end{aligned}
          :label: eqndefcombrule

Some or all combinations for different atom types can be given in the
``[ nonbond_params ]`` section, again with parameters V and
W as defined above. Any combination that is not given will be computed
from the parameters for the corresponding atom types, according to the
combination rule:

.. math:: \begin{aligned}
          \mbox{for combination rules 1 and 3}: & &
          \begin{array}{lll}
            C^{(6)}_{ij}  & = & \left(C^{(6)}_i\,C^{(6)}_j\right)^{\frac{1}{2}} \\
            C^{(12)}_{ij} & = & \left(C^{(12)}_i\,C^{(12)}_j\right)^{\frac{1}{2}}
          \end{array}
          \\
          \mbox{for combination rule 2}: & &
          \begin{array}{lll}
            \sigma_{ij}   & = & \frac{1}{2}(\sigma_i+\sigma_j) \\
            \epsilon_{ij} & = & \sqrt{\epsilon_i\,\epsilon_j}
          \end{array}\end{aligned}
          :label: eqngivencombrule

When :math:`\sigma` and :math:`\epsilon` need to be supplied (rules 2
and 3), it would seem it is impossible to have a non-zero :math:`C^{12}`
combined with a zero :math:`C^6` parameter. However, providing a
negative :math:`\sigma` will do exactly that, such that :math:`C^6` is
set to zero and :math:`C^{12}` is calculated normally. This situation
represents a special case in reading the value of :math:`\sigma`, and
nothing more.

There is only one set of combination rules for Buckingham potentials:

.. math:: \begin{array}{rcl}
          A_{ij}   &=& \left(A_{ii} \, A_{jj}\right)^{1/2}    \\
          B_{ij}   &=& 2 / \left(\frac{1}{B_{ii}} + \frac{1}{B_{jj}}\right)        \\
          C_{ij}   &=& \left(C_{ii} \, C_{jj}\right)^{1/2}
          \end{array}
          :label: eqnbuckinghamcombrule

Bonded parameters
~~~~~~~~~~~~~~~~~

The bonded
parameters
(*i.e.* bonds, bond angles, improper and proper dihedrals) are listed in
``ffbonded.itp``.  The entries in this database describe,
respectively, the atom types in the interactions, the type of the
interaction, and the parameters associated with that interaction. These
parameters are then read by
:ref:`grompp <gmx grompp>` when processing a
topology and applied to the relevant bonded parameters, *i.e.*
``bondtypes`` are applied to entries in the
``[ bonds ]`` directive, etc. Any bonded parameter that is
missing from the relevant :``[ *type ]`` directive generates
a fatal error. The types of interactions are listed in
:numref:`Table %s <tab-topfile2>`. Example excerpts from such files
follow:

::

    [ bondtypes ]
      ; i    j func        b0          kb
        C    O    1   0.12300     502080.
        C   OM    1   0.12500     418400.
        ......

    [ angletypes ]
      ; i    j    k func       th0         cth
       HO   OA    C    1   109.500     397.480
       HO   OA  CH1    1   109.500     397.480
       ......

    [ dihedraltypes ]
      ; i    l func        q0          cq
     NR5*  NR5    2     0.000     167.360
     NR5* NR5*    2     0.000     167.360
     ......

    [ dihedraltypes ]
      ; j    k func      phi0          cp   mult
        C   OA    1   180.000      16.736      2
        C    N    1   180.000      33.472      2
        ......

    [ dihedraltypes ]
    ;
    ; Ryckaert-Bellemans Dihedrals
    ;
    ; aj    ak      funct
    CP2     CP2     3       9.2789  12.156  -13.120 -3.0597 26.240  -31.495

In the ``ffbonded.itp`` file, you can add bonded parameters.
If you want to include parameters for new atom types, make sure you
define them in ``atomtypes.atp`` as well.

For most interaction types, bonded parameters are searched and assigned
using an exact match for all type names and allowing only a single set
of parameters. The exception to this rule are
dihedral
parameters. For
``[ dihedraltypes ]`` wildcard atom type names can be
specified with the letter ``X`` in one or more of the four
positions. Thus one can for example assign proper dihedral parameters
based on the types of the middle two atoms. The parameters for the entry
with the most exact matches, i.e. the least wildcard matches, will be
used. Note that |Gromacs| versions older than 5.1.3 used the first match,
which means that a full match would be ignored if it is preceded by an
entry that matches on wildcards. Thus it is suggested to put wildcard
entries at the end, in case someone might use a forcefield with older
versions of |Gromacs|. In addition there is a dihedral type 9 which adds
the possibility of assigning multiple dihedral potentials, useful for
combining terms with different multiplicities. The different dihedral
potential parameter sets should be on directly adjacent lines in the
``[ dihedraltypes ]`` section.
