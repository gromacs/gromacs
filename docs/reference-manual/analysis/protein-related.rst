Protein-related items
---------------------

| :ref:`gmx dssp <gmx dssp>`, :ref:`gmx rama <gmx rama>`,
  :ref:`gmx wheel <gmx wheel>`
| To analyze structural changes of a protein, you can calculate the
  radius of gyration or the minimum residue distances over time (see
  sec. :ref:`rg`), or calculate the RMSD (sec. :ref:`rmsd`).

You can also look at the changing of *secondary structure elements*
during your run. For this, you can use the program 
:ref:`gmx dssp <gmx dssp>`, which is a native implementation of DSSP algorithm :ref:`176 <refKabsch83>`. 

   Analysis of the secondary structure elements of a peptide in time.

One other important analysis of proteins is the so-called *Ramachandran
plot*. This is the projection of the structure on the two dihedral
angles :math:`\phi` and :math:`\psi` of the protein backbone, see
:numref:`Fig. %s <fig-phipsi>`: 

.. _fig-phipsi:

.. figure:: plots/phipsi.*
   :width: 8.00000cm

   Definition of the dihedral angles :math:`\phi` and :math:`\psi` of
   the protein backbone.

To evaluate this Ramachandran plot you can use the program
:ref:`gmx rama <gmx rama>`. A typical output
is given in :numref:`Fig. %s <fig-rama>`.

.. _fig-rama:

.. figure:: plots/rama.* 
    :width: 15.00000cm

    Ramachandran plot of a small protein.

When studying :math:`\alpha`-helices it is useful to have a *helical
wheel* projection of your peptide, to see whether a peptide is
amphipathic. This can be done using the :ref:`gmx wheel <gmx wheel>`
program. Two examples are plotted in
:numref:`Fig. %s <fig-hprwheel>`.

.. _fig-hprwheel:

.. figure:: plots/hpr-wheel.*
   :width: 10.00000cm

   Helical wheel projection of the N-terminal helix of HPr.

