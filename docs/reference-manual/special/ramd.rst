.. _ramd:

Random Acceleration Molecular Dynamics (RAMD)
---------------------------------------------

RAMD is a method to carry out molecular dynamics simulations with an additional randomly oriented
force applied to a molecule in the system. This is useful to accelerate the unbinding of ligands from
proteins or to study the egress pathways of ligands from proteins. The method is described in ref. \
:ref:`199 <refKokh2020>`.


Using RAMD
^^^^^^^^^^

RAMD simulations are enabled by the following :ref:`mdp` file options: ``ramd-active``.

Setting :mdp-value:`ramd-active = true` enables RAMD, using a configuration that can be
defined by specifying a RAMD configuration file using :mdp-value:`ramd-configfile`.

See :ref:`this section of the documentation <mdp-ramd>` for detailed usage of these options.


Configuration files for input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. mdp:: ramd-group1-force

   (600) [kJ mol\ :sup:`-1` nm\ :sup:`-1`]
   The force constant.
   For a set of compounds with the dissociation rate expected to vary
   within the range of 0.1-0.0001 1/s, a random force magnitude of 600
   kJ/mol/nm can be applied. If necessary, the force magnitude can be
   adjusted according to the longest and shortest dissociation time
   observed in simulations. The upper threshold of the force magnitude is
   determined by the fast-dissociated compounds, whose dissociation time
   should be longer than 100 ps. The lower threshold of the force magnitude
   depends on the computation facilities available.

.. mdp:: ramd-group1-r-min-dist

   (0.0025) [nm]
   This parameter affect absolute dissociation time but have less
   effect on the relative dissociation times of different compounds. It is
   recommended to use default value.

.. mdp:: ramd-group1-max-dist

   (4) [nm]
   This value has to be adjusted for the system studied: no
   protein-ligand contacts should be observed in the last snapshot of a
   dissociation trajectory. Usually 4 nm is enough, but in the case
   of a long dissociation channel (as in membrane proteins) maxDist must be
   increased accordingly. Method performance is not very sensitive to the
   upper limit of this parameter since motion of the free ligand due to the
   external force is very fast (i.e. the last part of the trajectory, where
   ligand does not interact with the protein, has a negligible contribution
   to the observed dissociation time.

.. mdp:: ramd-group1-receptor

   (Protein)

.. mdp:: ramd-group1-ligand

   (INH)

.. mdp:: ramd-group1-receptor-pbcatom

   (0)
   The value will be forwarded to the associated pull group of the receptor.

.. mdp:: ramd-group1-ligand-pbcatom

   (0)
   The value will be forwarded to the associated pull group of the ligand.

Example RAMD configuration file:

::
    ramd-group {
        receptor Protein
        ligand Ligand1
    }
    ramd-group {
        receptor Protein
        receptor-pbcatom 123
        ligand Ligand2
        ligand-pbcatom 456
        force 200
        max-dist 2.2
        r-min-dist 0.05
    }
