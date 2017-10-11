System preparation
==================

.. toctree::
   :hidden:

There are many ways to prepare a simulation system to run with
|Gromacs|. These often vary with the kind of scientific question being
considered, or the model physics involved. A protein-ligand atomistic
free-energy simulation might need a multi-state topology, while a
coarse-grained simulation might need to manage defaults that suit
systems with higher density.

Steps to consider
-----------------

The following general guidance should help with planning successful
simulations. Some stages are optional for some kinds of simulations.

1. Clearly identify the property or phenomena of interest to be
   studied by performing the simulation. Do not continue further until
   you are clear on this! Do not run your simulation and then seek to
   work out how to use it to test your hypothesis, because it may be
   unsuitable, or the required information was not saved.

2. Select the appropriate tools to be able to perform the simulation
   and observe the property or phenomena of interest. It is important
   to read and familiarize yourself with publications by other
   researchers on similar systems. Choices of tools include:

   - software with which to perform the simulation (consideration of
     force field may influence this decision)

   - the force field, which describes how the particles within the
     system interact with each other. Select one that is appropriate
     for the system being studied and the property or phenomena of
     interest. This is a very important and non-trivial step! Consider
     now how you will analyze your simulation data to make your
     observations.

3. Obtain or generate the initial coordinate file for each molecule to
   be placed within the system. Many different software packages are
   able to build molecular structures and assemble them into suitable
   configurations.

4. Generate the raw starting structure for the system by placing the
   molecules within the coordinate file as appropriate. Molecules may
   be specifically placed or arranged randomly. Several non-|Gromacs|
   tools are useful here; within |Gromacs| :ref:`gmx solvate`,
   :ref:`gmx insert-molecules` and :ref:`gmx genconf` solve frequent
   problems.

5. Obtain or generate the topology file for the system, using (for
   example) :ref:`gmx pdb2gmx`, :ref:`gmx x2top`, `SwissParam
   <http://swissparam.ch/>`_ (for CHARMM forcefield), `PRODRG
   <http://davapc1.bioch.dundee.ac.uk/cgi-bin/prodrg>`_ (for GROMOS96
   43A1), `Automated Topology Builder
   <http://compbio.biosci.uq.edu.au/atb/>`_ (for GROMOS96 53A6),
   `MKTOP <http://www.aribeiro.net.br/mktop>`_ (for OPLS/AA) or your
   favourite text editor in concert with chapter 5 of the |Gromacs|
   `Reference Manual`_. For the AMBER force fields, `antechamber
   <http://amber.scripps.edu/antechamber/antechamber.html>`__ or
   `acpype <https://github.com/choderalab/mmtools/blob/master/converters/acpype.py>`__
   might be appropriate.

6. Describe a simulation box (e.g. using :ref:`gmx editconf`) whose
   size is appropriate for the eventual density you would like, fill
   it with solvent (e.g. using :ref:`gmx solvate`), and add any
   counter-ions needed to neutralize the system (e.g. using :ref:`gmx
   grompp` and :ref:`gmx insert-molecules`). In these steps you may
   need to edit your topology file to stay current with your
   coordinate file.

7. Run an :ref:`energy minimization <gmx-energy-min>` 
   on the system (using :ref:`gmx grompp`
   and :ref:`gmx mdrun`). This is required to sort out any bad
   starting structures caused during generation of the system, which
   may cause the production simulation to crash. It may be necessary
   also to minimize your solute structure in vacuo before introducing
   solvent molecules (or your lipid bilayer or whatever else). You
   should consider using flexible water models and not using bond
   constraints or frozen groups. The use of position restraints and/or
   distance restraints should be evaluated carefully.

8. Select the appropriate simulation parameters for the equilibration
   simulation (defined in :ref:`mdp` file). You need to choose simulation
   parameters that are consistent with how force field was
   derived. You may need to simulate at NVT with position restraints
   on your solvent and/or solute to get the temperature almost right,
   then relax to NPT to fix the density (which should be done with
   Berendsen until after the density is stabilized, before a further
   switch to a barostat that produces the correct ensemble), then move
   further (if needed) to reach your production simulation ensemble
   (e.g. NVT, NVE). If you have problems here with the system :ref:`blowing
   up <blowing-up>`,
   consider using the suggestions on that page, e.g. position
   restraints on solutes, or not using bond constraints, or using
   smaller integration timesteps, or several gentler heating stage(s).

9. Run the equilibration simulation for sufficient time so that the
   system relaxes sufficiently in the target ensemble to allow the
   production run to be commenced (using :ref:`gmx grompp` and
   :ref:`gmx mdrun`, then :ref:`gmx energy` and `trajectory
   visualization tools
   <http://www.gromacs.org/Documentation/How-tos/Trajectory_Visualization>`_).

10. Select the appropriate simulation parameters for the production
    simulation (defined in :ref:`mdp` file). In particular, be careful not
    to re-generate the velocities. You still need to be consistent
    with how the force field was derived and how to measure the
    property or phenomena of interest.

.. _Reference Manual: `gmx-manual-parent-dir`_

Tips and tricks
---------------

Database files
^^^^^^^^^^^^^^

The ``share/top`` directory of a |Gromacs| installation contains
numerous plain-text helper files with the ``.dat`` file extension.
Some of the command-line tools (see :doc:`cmdline`) refer to these,
and each tool documents which files it uses, and how they are used.

If you need to modify these files (e.g. to introduce new atom types
with VDW radii into ``vdwradii.dat``), you can copy the file from your
installation directory into your working directory, and the |Gromacs|
tools will automatically load the copy from your working directory
rather than the standard one. To suppress all the standard
definitions, use an empty file in the working directory.
