.. _gmx-ff-included:

Force fields in |Gromacs|
=========================

.. _gmx-amber-ff:

AMBER
^^^^^

`AMBER`_ (Assisted Model Building and Energy Refinement) refers both to a set of molecular mechanical
:ref:`force fields <gmx-force-field>` for the simulation of biomolecules and a package of molecular simulation programs.

|Gromacs| supports the following AMBER force fields natively:

* AMBER94
* AMBER96
* AMBER99
* AMBER99SB
* AMBER99SB-ILDN
* AMBER03
* AMBERGS

Information concerning the force field can be found using the following information:

* `AMBER Force Fields <https://ambermd.org/AmberModels.php>`__ - background about the AMBER force fields
* `AMBER Programs <https://ambermd.org/AmberTools.php>`__ - information about the AMBER suite of
  programs for molecular simulation
* `ANTECHAMBER/GAFF <http://ambermd.org/antechamber/antechamber.html>`__ -
  Generalized Amber Force Field (GAFF) which is supposed to provide parameters
  suitable for small molecules that are compatible with the AMBER protein/nucleic
  acid force fields. It is available either together with AMBER, or through the
  antechamber package, which is also distributed separately. There are scripts
  available for converting AMBER systems (set up, for example, with GAFF) to
  |Gromacs| (`amb2gmx.pl <https://github.com/choderalab/mmtools/blob/master/converters/amb2gmx.pl>`__,
  or `ACPYPE <https://github.com/alanwilter/acpype>`_), but they do require 
  `AmberTools <https://ambermd.org/AmberTools.php>`_ installation to work.

.. _AMBER: http://ambermd.org/

.. _gmx-charmm-ff:

CHARMM
^^^^^^

`CHARMM`_ (Chemistry at HARvard Macromolecular Mechanics) is a both a set of force fields and
a software package for :ref:`molecular dynamics <gmx-md>` simulations and analysis. Includes united atom
(CHARMM19) and all atom (CHARMM22, CHARMM27, CHARMM36) :ref:`force fields <gmx-force-field>`.  The CHARMM27 force field
has been ported to |Gromacs| and is officially supported.  CHARMM36 force field files can be
obtained from the `MacKerell lab website`_, which regularly produces up-to-date CHARMM force field files in |Gromacs| format.

.. _CHARMM: http://www.charmm.org/
.. _MacKerell lab website: http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs

For using CHARMM36 in |Gromacs|, please use the following settings in the :ref:`mdp` file::

    constraints = h-bonds
    cutoff-scheme = Verlet
    vdwtype = cutoff
    vdw-modifier = force-switch
    rlist = 1.2
    rvdw = 1.2
    rvdw-switch = 1.0
    coulombtype = PME
    rcoulomb = 1.2
    DispCorr = no

Note that dispersion correction should be applied in the case of lipid monolayers, but not bilayers.

Please also note that the switching distance is a matter of some debate in lipid bilayer simulations,
and it is dependent to some extent on the nature of the lipid. Some studies have found that an 0.8-1.0 nm
switch is appropriate, others argue 0.8-1.2 nm is best, and yet others stand by 1.0-1.2 nm. The user
is cautioned to thoroughly investigate the force field literature for their chosen lipid(s) before beginning a simulation!

.. _gmx-gromos-ff:

GROMOS
^^^^^^

`GROMOS`_ is is a general-purpose molecular dynamics computer simulation package for the
study of biomolecular systems. It also incorporates its own force field covering proteins,
nucleotides, sugars etc. and can be applied to chemical and physical systems ranging from
glasses and liquid crystals, to polymers and crystals and solutions of biomolecules.

|Gromacs| supports the GROMOS force fields, with all parameters provided in the distribution
for 43a1, 43a2, 45a3, 53a5, 53a6 and 54a7. The GROMOS force fields are
:ref:`united atom force fields <gmx-force-field>`, i.e. without explicit aliphatic (non-polar) hydrogens.

* GROMOS 53a6 - in |Gromacs| format (J. Comput. Chem. 2004 vol. 25 (13): 1656-1676).
* GROMOS 53a5 - in |Gromacs| format (J. Comput. Chem. 2004 vol. 25 (13): 1656-1676).
* GROMOS 43a1p - 43a1 modified to contain SEP (phosphoserine), TPO (phosphothreonine),
  and PTR (phosphotyrosine) (all PO42- forms), and SEPH, TPOH, PTRH (PO4H- forms).

.. todo:: Add new force fields to the list

.. _GROMOS: https://www.igc.ethz.ch/gromos.html
.. _reference manual: gmx-manual-parent-dir_


.. _gmx-opls:

OPLS
^^^^

OPLS (Optimized Potential for Liquid Simulations) is a set of force fields developed by
Prof. William L. Jorgensen for condensed phase simulations, with the latest version
being `OPLS-AA/M <http://zarbi.chem.yale.edu/oplsaam.html>`__.

The standard implementations for those force fields are the *BOSS* and *MCPRO*
programs developed by the `Jorgensen group <http://zarbi.chem.yale.edu/software.html>`__

As there is no central web-page to point to, the user is advised to consult the
original literature for the `united atom (OPLS-UA) <https://doi.org/10.1021%2Fja00214a001>`__
and `all atom (OPLS-AA) <https://doi.org/10.1021%2Fja9621760>`__ force fields, as well as the
Jorgensen group `page <http://zarbi.chem.yale.edu/>`__
