Answers to frequently asked questions (FAQs)
============================================



.. Migrated from old website

.. toctree::
   :maxdepth: 2
   :hidden:

Questions regarding |Gromacs| installation
------------------------------------------

1. Do I need to compile all utilities with parallelization?
   
   In general only the :ref:`mdrun <gmx mdrun>` binary is able to use the parallelization 
   offered by :ref:`MPI <mpi-support>`. So you only need to use the ``-DGMX_MPI=on`` flag
   when :ref:`configuring <configure-cmake>` in combination with
   ``-DGMX_BUILD_MDRUN_ONLY=ON``.


2. Should my version be compiled using double precision?

   In general |Gromacs| only needs to be build in single precision, as random effects
   due to :ref:`floating point arithmetic <gmx-floating-point>` will generate random
   behaviour after a certain time regardless of precision. Still, you should consider your needs
   when using |Gromacs|, and can use double precision builds for e.g.
   :ref:`accurate minimization <gmx-energy-min>`. Other usage my also depend on your target system
   and should be decided upon according to the :ref:`individual instructions <gmx-special-build>`.

Questions concerning system preparation and preprocessing
---------------------------------------------------------

1. Where can I find a solvent :ref:`coordinate file <gmx-structure-files>` for use with :ref:`solvate <gmx solvate>`?

   The files defining solvent types to be used with :ref:`solvate <gmx solvate>` are found in the respective trees
   of the :ref:`force field <gmx-force-field>`, with the :ref:`structure files <gmx-structure-files>` found 
   in the ``$GMXDIR/share/gromacs/top`` directory.    Those files contain equilibrated boxes of solvent that
   are recognized directly by |Gromacs| when using e.g. ``-cs spc216.gro``    as the argument to
   :ref:`solvate <gmx solvate>`. Other solvent boxes can be prepared by the user as described in
   :ref:`Non-Water Solvation <gmx-solvate-other>` and on the manual page for :ref:`solvate <gmx solvate>`.

2. How to prevent :ref:`solvate <gmx solvate>` from placing waters in undesired places?

   Water placement is in general well behaved when solvating proteins, but can be difficult when setting up
   :ref:`membrane <gmx-membrane>` or micelle simulations. In those cases, waters my be placed in between the
   alkyl chains of the lipids, leading to problems later :ref:`during the simulation <blowing-up>`.
   You can either remove those waters by hand (and do the accounting in the :ref:`topology <top>` file),
   or set up a local, membrane simulation specific copy of the ``vdwradii.dat`` file from the ``$GMXLIB``
   directory for your project, with an increased :ref:`vdW radius<gmx-vdw>`. Recommended e.g.
   at a common `tutorial`_ is the use of 0.375 instead of 0.15.

.. _tutorial: http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/membrane_protein/03_solvate.html

3. How do I provide multiple definitions of bonds / dihedrals in a topology?

   You can add additional bonded terms beyond those that are normally defined for a residue (e.g. when defining
   a special ligand) by including additional copies of the respective lines under the
   ``[ bonds ]``, ``[ pairs ]``, ``[ angles ]`` and ``[ dihedrals ]`` sections in either the :ref:`itp` file
   or the final :ref:`topology <top>`. This will **add** those extra terms to the potential energy evaluation,
   but **will not** remove the previous ones. So be careful with duplicate entries. Also keep in mind that this **does not** mean
   duplicated entries for ``[ bondtypes ]``, ``[ angletypes ]``, or ``[ dihedraltypes ]``, where duplicates overwrite the
   previous values.

.. old page says this is not the case when using CHARMM, is this still the case?


4. Do I really need a :ref:`gro` file?

   The :ref:`gro` file is used in |Gromacs| as a unified :ref:`structure file <gmx-structure-files>` format
   that can be read by all processes. The large majority of |Gromacs| routines can also use other file
   types such as :ref:`pdb`, with the limitations that no velocities are available in :ref:`this case <gmx-need-for-gro>`.

5. Do I always need to run :ref:`pdb2gmx <gmx pdb2gmx>` when I already produced an :ref:`itp` file elsewhere (like `PRODRG`_)?

   You don't need to prepare additional files if you already have all :ref:`itp` and :ref:`top` files prepared through other tools.

.. _PRODRG: http://davapc1.bioch.dundee.ac.uk/cgi-bin/prodrg

6. How can I build in missing atoms?

   |Gromacs| has no internal methods on how to add missing atoms (except hydrogens). If your system is missing some part,
   you will have to add the missing pieces using external programs to avoid the :ref:`missing atom <gmx-atom-missing>`
   error. This can be done using programs such as `Chimera <https://www.cgl.ucsf.edu/chimera/>`__ in combination
   with `Modeller <https://salilab.org/modeller/>`__, `Swiss PDB Viewer <https://spdbv.vital-it.ch/>`__,
   `Maestro <https://www.schrodinger.com/maestro>`__. **Do not run** a simulation that had missing atoms except
   if you know exactly what you are doing.

7. Why is the total charge of my system not an integer like it should be?

   In :ref:`floating point <gmx-floating-point>` math, real numbers can not be displayed to arbitrary precision
   (for more on this see e.g. `Wikipedia <https://en.wikipedia.org/wiki/Floating-point_arithmetic>`__). This means
   that very small differences to the final integer value will persist, and |Gromacs| will not lie to you and
   round those values up or down. If your charge differs from the integer value by a larger amount, this usually means that
   something went wrong during your system preparation
