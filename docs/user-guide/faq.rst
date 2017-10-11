Answers to frequently asked questions (FAQs)
============================================

.. _reference manual: `gmx-manual-parent-dir`_

.. Migrated from old website

.. toctree::
   :maxdepth: 2
   :hidden:

Questions regarding |Gromacs| installation
------------------------------------------

#. Do I need to compile all utilities with parallelization?
   
   In general only the :ref:`mdrun <gmx mdrun>` binary is able to use the parallelization 
   offered by :ref:`MPI <mpi-support>`. So you only need to use the ``-DGMX_MPI=on`` flag
   when :ref:`configuring <configure-cmake>` in combination with
   ``-DGMX_BUILD_MDRUN_ONLY=ON``.


#. Should my version be compiled using double precision?

   In general |Gromacs| only needs to be build in single precision, as random effects
   due to :ref:`floating point arithmetic <gmx-floating-point>` will generate random
   behaviour after a certain time regardless of precision. Still, you should consider your needs
   when using |Gromacs|, and can use double precision builds for e.g.
   :ref:`accurate minimization <gmx-energy-min>`. Other usage my also depend on your target system
   and should be decided upon according to the :ref:`individual instructions <gmx-special-build>`.

Questions concerning system preparation and preprocessing
---------------------------------------------------------

#. Where can I find a solvent :ref:`coordinate file <gmx-structure-files>` for use with :ref:`solvate <gmx solvate>`?

   The files defining solvent types to be used with :ref:`solvate <gmx solvate>` are found in the respective trees
   of the :ref:`force field <gmx-force-field>`, with the :ref:`structure files <gmx-structure-files>` found 
   in the ``$GMXDIR/share/gromacs/top`` directory.    Those files contain equilibrated boxes of solvent that
   are recognized directly by |Gromacs| when using e.g. ``-cs spc216.gro``    as the argument to
   :ref:`solvate <gmx solvate>`. Other solvent boxes can be prepared by the user as described in
   :ref:`Non-Water Solvation <gmx-solvate-other>` and on the manual page for :ref:`solvate <gmx solvate>`.

#. How to prevent :ref:`solvate <gmx solvate>` from placing waters in undesired places?

   Water placement is in general well behaved when solvating proteins, but can be difficult when setting up
   :ref:`membrane <gmx-membrane>` or micelle simulations. In those cases, waters my be placed in between the
   alkyl chains of the lipids, leading to problems later :ref:`during the simulation <blowing-up>`.
   You can either remove those waters by hand (and do the accounting in the :ref:`topology <top>` file),
   or set up a local, membrane simulation specific copy of the ``vdwradii.dat`` file from the ``$GMXLIB``
   directory for your project, with an increased :ref:`vdW radius<gmx-vdw>`. Recommended e.g.
   at a common `tutorial`_ is the use of 0.375 instead of 0.15.

.. _tutorial: http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/membrane_protein/03_solvate.html

#. How do I provide multiple definitions of bonds / dihedrals in a topology?

   You can add additional bonded terms beyond those that are normally defined for a residue (e.g. when defining
   a special ligand) by including additional copies of the respective lines under the
   ``[ bonds ]``, ``[ pairs ]``, ``[ angles ]`` and ``[ dihedrals ]`` sections in either the :ref:`itp` file
   or the final :ref:`topology <top>`. This will **add** those extra terms to the potential energy evaluation,
   but **will not** remove the previous ones. So be careful with duplicate entries. Also keep in mind that this **does not** mean
   duplicated entries for ``[ bondtypes ]``, ``[ angletypes ]``, or ``[ dihedraltypes ]``, where duplicates overwrite the
   previous values.

.. old page says this is not the case when using CHARMM, is this still the case?


#. Do I really need a :ref:`gro` file?

   The :ref:`gro` file is used in |Gromacs| as a unified :ref:`structure file <gmx-structure-files>` format
   that can be read by all processes. The large majority of |Gromacs| routines can also use other file
   types such as :ref:`pdb`, with the limitations that no velocities are available in :ref:`this case <gmx-need-for-gro>`.

#. Do I always need to run :ref:`pdb2gmx <gmx pdb2gmx>` when I already produced an :ref:`itp` file elsewhere (like `PRODRG`_)?

   You don't need to prepare additional files if you already have all :ref:`itp` and :ref:`top` files prepared through other tools.

   Examples for those are `CHARMM-GUI <http://www.charmm-gui.org/>`__, `ATB (Automated Topology Builder <https://atb.uq.edu.au/>`__
   and `pmx <http://pmx.mpibpc.mpg.de/instructions.html>`__.

.. _PRODRG: http://davapc1.bioch.dundee.ac.uk/cgi-bin/prodrg

#. How can I build in missing atoms?

   |Gromacs| has no internal methods on how to add missing atoms (except hydrogens). If your system is missing some part,
   you will have to add the missing pieces using external programs to avoid the :ref:`missing atom <gmx-atom-missing>`
   error. This can be done using programs such as `Chimera <https://www.cgl.ucsf.edu/chimera/>`__ in combination
   with `Modeller <https://salilab.org/modeller/>`__, `Swiss PDB Viewer <https://spdbv.vital-it.ch/>`__,
   `Maestro <https://www.schrodinger.com/maestro>`__. **Do not run** a simulation that had missing atoms except
   if you know exactly what you are doing.

#. Why is the total charge of my system not an integer like it should be?

   In :ref:`floating point <gmx-floating-point>` math, real numbers can not be displayed to arbitrary precision
   (for more on this see e.g. `Wikipedia <https://en.wikipedia.org/wiki/Floating-point_arithmetic>`__). This means
   that very small differences to the final integer value will persist, and |Gromacs| will not lie to you and
   round those values up or down. If your charge differs from the integer value by a larger amount, this usually means that
   something went wrong during your system preparation

Questions regarding simulation methodology
------------------------------------------

#.  Should I temperature-couple a handful of ions to their own thermal bath?

    Short answer: **No**. Long answer: You need to consider the minimal size of your
    temperature coupling groups, as explained in :ref:`gmx-thermostats` and more
    specifically in :ref:`gmx-thermostats-dont`, as well as the implementation
    of your chosen thermostat as described in the `reference manual`_.

#.  Why do my grompp restarts always start from time zero?

    You need to set the correct options in the :ref:`mdp` file to ensure the
    correct continuation, as explained in the `reference manual`_.

.. TODO make links work :ref:`Continuing simulations <gmx-cont-simulation>`.

#.  Why can't I do conjugate gradient minimization with constraints?

    Minimization with the conjugate gradient scheme can not be performed with constrains
    as described in the `reference manual`_, and some additional information
    on `Wikipedia <https://en.wikipedia.org/wiki/Conjugate_gradient_method>`__.

#.  Why can't I constrain atoms in place in my energy minimization or simulation?

    In general a ``constrain`` would mean that the atoms are fixed in space, something 
    that is not possible in |Gromacs|. What is possible is using a set of position
    restraints to place movement penalties on the atoms. Files that control this
    kind of behaviour can be created using :ref:`genrestr <gmx genrestr>`.

#.  How do I extend a completed a simulation to longer times?

    Please see the section on Continuing simulations in the `reference manual`_.
    You can either prepare a new :ref:`mdp` file or extend the simulation time
    in the original :ref:`tpr` file using :ref:`convert-tpr`<gmx convert-tpr>`.

.. TODO #.  How do I complete a crashed simulation?

..    This can be easily achieved using the checkpoint reading
    :ref:`available <gmx-cont-crash>` in |Gromacs| versions newer than 4.

.. TODO #.  How can I do a simulation at constant pH?

..    This is a rather large topic, and you should at least read the short
    :ref:`Constant pH How-To <gmx-howto-cph>` and all of the literature
    included there to get an overview over the topic.

.. TODO #.  How should I compute a single-point energy?

..    This can be achieved with the ``-rerun`` option to :ref:`mdrun <gmx mdrun>`.
    Please see the description of the command line help and the
    :ref:`single point <gmx-single-point-energy>` explanation for further information.

Parameterization and Force Fields
---------------------------------

#.  I want to simulate a molecule (protein, DNA, etc.) which complexes with
    various transition metal ions, iron-sulfur clusters, or other exotic species.
    Parameters for these exotic species aren't available in force field X.
    What should I do?

    First, you should consider on how well :ref:`MD <gmx-md>` will actually describe your
    system (e.g. see some of the `recent literature <https://dx.doi.org/10.1021%2Facs.chemrev.6b00440>`__).
    Then you need to prepare your own set of parameters and add a new residue
    to your :ref:`force field <gmx-force-field>` of choice. Then you will have to validate that
    your system behaves in a physical way, before continuing your studies. You could
    also try to build a more simplified model that does not rely on the complicated additions,
    as long as it still represents the correct *real* object in the laboratory.

#.  Should I take parameters from one force field and apply them inside another that is missing them?

    Short answer: **NO**. Long answer: No! Molecules parametrized for a given
    :ref:`force field <gmx-force-field>` will not behave in a physical manner when interacting with
    other molecules that have been parametrized according to different standards. If your
    required molecule is not included in the force field you need to use, you will 
    have to parametrize it yourself according to the methodology of this force field.

Analysis and Visualization
--------------------------

.. TODO #.  How do I visualize a trajectory?

..    Use one of the number of different programs that can visualize
    coordinate :ref:`files and trajectories <gmx-howto-visualize>`.

#.  Why am I seeing bonds being created when I watch the trajectory?

    Most visualization software only analyse the bond status of atoms depending
    on a set of predefined distances. So the bonding pattern created by them
    might not be the one defined in your :ref:`topology <top>` file. What
    matters is the information encoded in there.

#.  When visualizing a trajectory from a simulation using PBC, why are there holes or my peptide leaving the simulation box?

    Those holes and molecules moving around are just a result of molecules
    ranging over the :ref:`box boundaries and wrap around <gmx-pbc>`
    and are not a reason for concern. You can fix the visualization using :ref:`trjconv <gmx trjconv>`
    to prepare the structure for analysis.

#.  Why is my total simulation time not an integer like it should be?

    As the simulation time is calculated using :ref:`floating point arithmetic <gmx-floating-point>`,
    rounding errors can occurs but are not of concern.


