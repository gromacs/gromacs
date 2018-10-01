.. _constraintalg:

Constraint algorithms
---------------------

Constraints are defined in the ``[ constraints ]`` section. The format is two atom numbers
followed by the function type, which can be 1 or 2, and the constraint
distance. The only difference between the two types is that type 1 is
used for generating exclusions and type 2 is not (see sec. :ref:`excl`).
The distances are constrained using the LINCS or the SHAKE algorithm,
which can be selected in the :ref:`mdp` file. Both types of constraints can be
perturbed in free-energy calculations by adding a second constraint
distance (see :ref:`constraintforce`). Several types of bonds and
angles (see :numref:`Table %s <tab-topfile2>`) can be converted automatically to
constraints by :ref:`grompp <gmx grompp>`. There are several options for this in the :ref:`mdp`
file.

We have also implemented the SETTLE
algorithm \ :ref:`47 <refMiyamoto92>`, which is an analytical solution of SHAKE, specifically for
water. SETTLE can be selected in the topology file. See, for instance,
the SPC molecule definition:

::

    [ moleculetype ]
    ; molname       nrexcl
    SOL             1

    [ atoms ]
    ; nr    at type res nr  ren nm  at nm   cg nr   charge
    1       OW      1       SOL     OW1     1       -0.82
    2       HW      1       SOL     HW2     1        0.41
    3       HW      1       SOL     HW3     1        0.41

    [ settles ]
    ; OW    funct   doh     dhh
    1       1       0.1     0.16333

    [ exclusions ]
    1       2       3
    2       1       3
    3       1       2

The ``[ settles ]`` directive defines the first atom of the
water molecule. The settle funct is always 1, and the distance between
O-H and H-H distances must be given. **Note** that the algorithm can
also be used for TIP3P and TIP4P \ :ref:`128 <refJorgensen83>`. TIP3P just has
another geometry. TIP4P has a virtual site, but since that is generated
it does not need to be shaken (nor stirred).
