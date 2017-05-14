=============
Python module
=============

Goal: Expose Gromacs functionality in an interface
that can be easily driven and interacted with from Python (or other
high-level languages) and external APIs for e.g. "cloud" computing resources,
container systems, work flow / compute resource managers. Data graph execution
at various levels may be deferred to other tools, but should at least be
facilitated by the design.

There are several broad classes of use cases that warrant different levels of access,
abstraction, performance, and invasiveness into the existing code.
In rough order of increasing invasiveness.

1. Post-run analysis (expose and extend the current C++ API, file readers, wrapped tools)
2. Preparation, preprocessing, and editing of simulation jobs or on-disk information (encapsulation of tools like `pdb2gmx`, `grompp`, `trajedit`, `solvate`, and of data structures such as topologies, `tpr` files, trajectories and checkpoints)
3. Controlling Gromacs code (scripting, managing compute resources, high level workflow and data flow, clearly separating model parameters / simulation parameters / workflow parameters / runtime parameters)
4. Live analysis (hooks into MD loop: callbacks, custom loggers)
5. Extension of Gromacs code (adding custom potentials, logic, low level workflow and data flow)

The APIs developed should allow access to these functionalities.
A user should be able to link together tasks programatically, dynamically, and
with no requirements of filesystem I/O.


Gromacs command-line functionality should be easily achieved in a Python script. This
does not mean that all CLI tools should be simply exposed to Python, since
many represent trivial data transformations, tasks that are well-handled by other
Python-based tools, or easily-scripted tasks. Common scriptable tasks can be
provided by high-level convenience functions or "bottled" workflows. More
elemental functionality that should be available through the C++ API can be
exposed through bindings to the C++ API.

With Python's reference counting and conceptual differences about the mutability
of data, it is appropriate to support a few different modes of access to Gromacs
data structures. When it is necessary to copy data, explicit API syntax should
make it clear to the user. When possible, the Python buffer protocol can be used
to give Python direct access to memory shared by C++ data structures. In general,
though, for data and objects that exist only to be used in subsequent API calls,
proxy objects can be used at a high level with interactions between objects
deferred to the C++ library API where the parallelism framework and optimizations
already exist. This also helps with the notion that the Python interpreter
expects to "own" objects that are handed to it through bindings.

The design should emphasize idiomatic Python usage.

Suggested use case timeline:

1. Read and write various Gromacs files with appropriately
general objects at the high level interface.
2. Integrate with the Trajectory Analysis framework to run pipelines of modules
implemented in C++, Python, and external tools.
3. Drive a simulation workflow from a Python script.
4. Manage an adaptive simulation workflow combining analysis and simulations
configured programmatically, with minimal I/O or parallel gathers/scatters and
with workflow checkpointing. Edit input configurations, topologies, and runtime parameters.
5. Use Python call-backs for low-performance invasive experimentation with a running
simulation and provide the user interface to connect plugins that use the C++
plugin API.

Other Python Gromacs interfaces
===============================

`gmxscript <https://github.com/pslacerda/gmx>`__ : all-python CLI
wrapper manages files, provides "checkpointed" re-runnable workflow with
mdp rewriting and CLI replacement

`GromacsWrapper <http://gromacswrapper.readthedocs.io/en/latest/>`__
(Beckstein) : all-python CLI wrapper provides thorough scriptable
interface with error handling

`pmx <https://github.com/dseeliger/pmx>`__ (D. Seeliger) : Manipulate
Gromacs files and data structures

`grompy <https://github.com/GromPy>`__ (Rene Pool) : patches old Gromacs
code to provide ctypes interface to Gromacs libraries

`gmx\_wrapper <https://github.com/khuston/gmx_wrapper>`__ : very
bare-bones CLI wrapper

`GromacsPipeline <https://biod.pnpi.spb.ru/gitweb/?p=alexxy/gromacs.git;a=commit;h=1241cd15da38bf7afd65a924100730b04e430475>`__
(Redmine 1625) adds SIP bindings to run Gromacs analysis modules in a
single pass on a trajectory

Existing Python tools to leverage
=================================

| Some of the tools available from the ``gmx`` command-line interface
  are bundled largely for convenience,
| and Python users may be better served by turning to other projects
  rather than relying on exposing
| redundant functionality from ``libgromacs``.
| In other cases, Gromacs developers may simply want to be aware of the
  other tools and interoperability requirements.

-  `pypdb <https://github.com/williamgilpin/pypdb>`__ A Python API for
   the RCSB Protein Data Bank (PDB)
-  `chemlab <http://chemlab.github.io/chemlab/>`__ Qt-powered molecular
   viewer
-  `chemview <https://github.com/gabrielelanaro/chemview>`__ molecular
   viewer for IPython from the author of chemlab
-  `MDTraj <http://mdtraj.org/>`__ Read, write and analyze MD
   trajectories
-  `MDAnalysis <http://www.mdanalysis.org>`__ analyze trajectories from
   molecular dynamics
-  `Biopython <https://github.com/biopython/biopython>`__ tools for
   computational molecular biology
-  several packages named "Biopy" do not seem relevant
-  `PyMOL <http://www.pymol.org/>`__ molecular viewer
-  `mmLib <http://pymmlib.sourceforge.net/>`__ Python Macromolecular
   Library analysis, manipulation, and viewing
-  `msmbuilder <http://msmbuilder.org/>`__ Statistical models for
   Biomolecular Dynamics
-  `PyEMMA <http://emma-project.org/>`__ estimation, validation and
   analysis Markov models from molecular dynamics (MD) data
-  `ExSTACY <http://extasy-project.org>`__ Collection of tools that
   couples large ensemble Molecular Dynamics calculations with analysis
   to provide improvements in sampling. Includes
-  `EnsembleMD <https://github.com/radical-cybertools/radical.ensemblemd>`__
   framework for running ensemble workflows
-  `COCO <https://bitbucket.org/extasy-project/coco>`__
-  `LSDMap <https://sourceforge.net/projects/lsdmap/>`__ Locally Scaled
   Diffusion Maps
-  `pyPcazip <https://bitbucket.org/ramonbsc/pypcazip>`__ principle
   component analysis
-  `Plumed <http://www.plumed.org>`__ free energy calculation using
   collective variables
-  pdb tools
-  OpenMM
-  msmbuilder
-  BioXL

Naming
======

| The package should have a name that does not collide with other known
  projects, particularly any projects on pypi.
| "gmx" is available on pypi, but similar existing package names include

-  ``gmx-script``
-  ``python-gmx``

``pygmx`` is the name used for a bindings package posted to the Gromacs
list that resides at https://biod. pnpi.spb.ru/~alexxy/pygmx

There are several packages named ``grompy``

The Beckstein lab has ``GromacsWrapper`` which exists on
PythonHosted.org, GitHub, and PyPi.

Python version
==============

Python 3.3 has a lot of improvements and changes, including better 2.7 compatibility

Python 3.4

* built-in pip
* enum types
* new pickle protocol

Python 3.5

* typing and coroutines
* RecursionError exception
* Generators have gi_yieldfrom
* memoryview tuple indexing
* hexadecimal methods

Linux distributions released after June 2013 and supported at least to June 2019.

| Ubuntu 14.04 (trusty): Python 3.4
| Ubuntu 16.04 (denial): Python 3.5
| Debian 8 (jessie): Python 3.4
| Debian 9 (stretch): Python 3.5
| Linux Mint 18 (rosa): 3.4
| Linux Mint 17 (sarah): 3.5
| Fedora 23: 3.4
| Fedora 24+: 3.5
| RHEL 7: n/a
| CentOS 7: n/a

Suggestion: require Python 3.4+: widely supported by distributions that, as of a
projected public release, will have been released less than five years ago and
will be supported for at least another year. Built-in `pip` and enum types would
be very nice to be able to rely on. If Python 2.7 is supported, then we should
not support less than Python 3.3, but it may be hard enough to take advantage of
features in 3.4+ that there isn't a good reason to require it without motivation
from a clear roadmap of which versions Gromacs will support in the future.

Use cases and Scenarios
=======================

Traditional HPC scenarios
-------------------------
An incomplete list of some contemporary modes of scientific computing.

large-scale MPI job on HPC cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Traditional large-scale MPI job on HPC cluster

|  1. User submits job requesting many resources (e.g. N m-core ranks)
|  2. Queuing system simultaneously launches N processes with access to MPI
    environments and interconnect libraries installed on the compute nodes
|  3. Each process uses its linked MPI library and runtime environment to access its MPI communicator
|  4. Each process runs a copy of the same object code with different(ly processed) input
|  5. Processes communicate through the MPI communicator as needed
|  6. Processes terminate individually (but hopefully simultaneously) and job
    completes, leaving results in a shared filesystem or first transferring data to a master rank

user-managed bag-of-jobs on HPC cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Traditional independent bag-of-jobs on HPC cluster with shared storage

|  1. User submits many isolated job scripts
|  2. Queuing system executes job script on reserved resources as available.
|  3. Job script uses unique parameters, environment variables, and job ID to run a unique simulation asynchronously.
|  4. User retrieves and/or analyzes data

Alternate scenarios

|  1b. User submits a job script many times (or as a job array) with different parameters
|  4b. Queuing system job dependencies chain tasks (unreliable in practice) using shared filesystem

pilot-managed jobs on HPC cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Traditional independent bag-of-jobs on pilot-managed HPC resource allocation

|  1. User submits pilot job
|  2. Queuing system notifies user when pilot job is running
|  3a. User submits bag-of-jobs to pilot manager
|  4a. Pilot manager transfers scripts and data and connects to compute nodes to execute tasks as resources are available in the pilot job

Alternate scenario: Managed ensemble workflow on pilot-managed HPC resource allocation

|  3b. User specifies workflow elements and parameters to ensemble manager
    (e.g. cross-correlate all permutations of trajectory pairs)
|  4b. Ensemble manager generates pilot work units using pilot interface
|  5. Pilot manager finds available resources, transfers data and connects to
    compute nodes to execute remote commands (which may invoke single processes,
    multi-threaded processes, or MPI process groups).
|  6. Ensemble manager collates results

commodity cloud computing
~~~~~~~~~~~~~~~~~~~~~~~~~
Various scenarious, depending on environment and user preferences, involving VMs
or Docker images. Avoiding expensive data storage costs can require a lot of
scripting and testing of chained tasks.

Target scenarios
----------------

Following are an assortment of tasks that we are targeting with the Python
interface and API design.

Programmatically launch many similar simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run -> analyze -> run chained tasks in pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data graph execution
~~~~~~~~~~~~~~~~~~~~

Implement/adapt a new execution manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Scripted workflow
~~~~~~~~~~~~~~~~~

| 1. `Load pdb data`_
| 2. `Perform energy minimization`_
| 3. `Solvate`_
| 4. `Add ions`_
| 5. `Perform energy minimization`_
| 6. `Equilibration`_
|   6.1. NVT
|   6.2. NpT
| 7. Launch simulation job

Load pdb data
~~~~~~~~~~~~~

| 1. Read and clean pdb file
| 2. Get topology from factory using built-in force field definition
| 3. Produce ``.gro`` and ``.top`` files
| 4. `Set up a simulation box`_

Alternate: bottled workflow

| 1a. use a utility / helper function ``like pdb2gmx``
|   1.1. proceed to 4.

Alternate: non-automated topology construction

| 2a. Start from an empty or existing Topology object
|   2a.1. get or define a forcefield data structure
|   2a.2. use member functions of Topology object to modify

Set up a simulation box
~~~~~~~~~~~~~~~~~~~~~~~

| 1. Read a structure (``.gro``) file
| 2. Edit data
| 3. Write out new structure file

Alternate: no file I/O

| 1a. Create or reuse a AtomData object
|   1a.1 set box parameters for object
|   1a.2 done

Perform energy minimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

| 1. Create input record
|   1.1 Read mdp file
|   1.2 Save tpr file
| 2. Optionally configure filesystem output and then run.

Alternate: avoid redundant file I/O

| 1a. Operate on and produce Topology, Particles, and System objects without
   file I/O using previously generated objects and mdp file
| 2a. Run energy minimizer with no filesystem output

Alternate: skip mdp file

| 1b. Skip mdp file with granular functionality and parameters in Python data structures
|   1.1 set parameters for minimizer, general integrator, and neighborlist
|   1.2 create System and Context objects
|   1.3 create minimizer object and bind everything together

or

| 1c. reuse the last minimizer System and Context
|   1c.1. update minimizer parameters and configuration
|   1c.2. proceed to 2.

Alternate: additional logging

| 2d. Optionally configure Reporter objects before running

Solvate
~~~~~~~

1. Load solvent and coordinates files
2. Use utility to create new configuration from appropriate data
3. Save solvated ``.gro`` file

Alternate: use package data instead of files

| 1a. Specify solvent molecule accessed through the module

Alternates: System methods instead of utility functions

| 2a. Update an existing System with provided solvent data
| 2b. Update an existing System with the default solvent

Add ions
~~~~~~~~

| 1. Get configuration with force field
|  1.1 Read new mdp file
|  1.2 Save new tpr file
| 2. Use utility to insert ions and produce a new System object
| 3. Optionally save new configuration

Alternate:

| 1a. re-use an instantiated System

User alternative:

| 2a. new system constructed with fancy constructor

Implementation alternative:

| 2b. utility produces just a configuration

Alternate: operate on objects

| 2c. Use generic ``change_residues`` method with various modes of operation.

Equilibration
~~~~~~~~~~~~~

| 1. Configure from mdp and other files using convenience function to
   construct System
| 2. Optionally save input record
| 3. Run simulation
| 4. Optionally save output to file or just retain object for later API call

Implementation option:

| 1a. convenience function can handle a variety of input argument types and
   forward to an appropriate helper function for files, current instances, etc.

Examples
==============

The following examples show what the Python interface might look like and
illustrate more concretely the concepts in the scenarios above.
Ultimately this documentation-by-example will be extracted from the source code,
but right now it is static content that will guide the source code layout.

Sea-level Scenario: Simulate counter-ion effects on peptide from PDB
--------------------------------------------------------------------

Reimplement the CLI workflow described in the funnel web spider toxin
tutorial at
http://cinjweb.umdnj.edu/~kerrigje/pdf_files/fwspidr_tutor.pdf

1. Prepare a configuration from PDB data
2. In vacuo energy minimization
3. Solvate
4. Add ions
5. Solvated energy minimization
6. Two-step equilibration

Assumes

.. sourcecode :: python

     import gmx

Prepare a configuration from PDB data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read and clean pdb file.

Consume ``fws.pdb`` and produces ``.gro`` and ``.top`` files. e.g.
``pdb2gmx -ignh -ff amber99sb-ildn -f fws.pdb -o fws.gro -p fws.top -water tip3p``

.. sourcecode :: python

     # Depending on how pdb2gmx is currently implemented...
     cleanpdb = gmx.util.clean_pdb('fws.pdb', ignore_hydrogens=True)
     (atoms, topology) = gmx.util.pdb2gmx(cleanpdb,
                             forcefield=gmx.forcefield.amber99sb-ildn)
     # or
     (atoms, topology) = gmx.util.pdb2gmx('fws.pdb',
                             forcefield=gmx.forcefield.amber99sb-ildn,
                             ignore_hydrogens=True)

     # Alternatively
     import pypdb
     pdbfile = pypdb.get_pdb_file('1OMB', filetype='pdb', compression=True)
     # Then use other tools to "clean" the PDB record

Build topology using force field and structure data

.. sourcecode :: python

     # return a gmx.Topology object for atoms in pdbfile
     topology = gmx.TopologyBuilder(gmx.forcefields.amber99sb-ildn).from_pdb(pdbfile)
     # then extend topology information with chosen water model
     topology.add_moleculetype('SOL', gmx.amber99sb-ildn.tip3p_flexible)

Produce .gro and .top files
with a convenience function.

.. sourcecode :: python

     gmx.AtomData.from_pdb(pdbfile).save('fws.gro')
     topology.save('fws.top')

Alternatively

.. sourcecode :: python

     # Use utility to read and write files. Bottled workflow.
     gmx.utils.pdb2gmx(force_field=gmx.forcefield.amber99sb-ildn,
                                     pdb_file="fws.pdb",
                                     coords_file="fws.gro",
                                     topology_file="fws.top",
                                     water=gmx.forcefield.amber99sb-ildn.tip3p)

Includes: Set up the box

Set up the box.
~~~~~~~~~~~~~~~

e.g. ``editconf -f fws.gro -o fws-PBC.gro -bt dodecahedron -d 1.2``

.. sourcecode :: python

     #Read the file into a Atom data object, intuiting the file type by extension
     grofile = gmx.AtomData('fws.gro')
     # Edit the atom data object
     grofile.set_box(gmx.SimulationBox(pbc='dodecahedron', L=1.2))
     # Write out file using the appropriate writer for the file extension
     grofile.save('fws-PBC.gro')

     # or

     # Edit the atom data object
     atoms.set_box(gmx.SimulationBox(pbc='dodecahedron', L=1.2))

In vacuo energy minimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prepare input record. e.g. ``grompp -f em-vac-pme.mdp -c fws-PBC.gro -p fws.top -o em-vac.tpr``

.. sourcecode :: python

     # using convenience function for mdp file
     minimization = gmx.System.from_mdpfile(mdpfile='em-vac-pmd.mdp',
                                                     topology=topology,
                                                     atoms=grofile,
                                                     )
     # or
     minimization = gmx.System.from_mdpfile(mdpfile='em-vac-pmd.mdp',
                                                     topology=topology,
                                                     atoms=atoms,
                                                     )
     # optionally
     # save initialization file
     minimization.save('em-vac.tpr')

More granularly, we might use key-value parameters and access objects more directly.
periodic boundary conditions are already defined in the configuration box
no bond types are replaced by constraints, so we could ignore "constraints=none".
If we were to replace bond types with constraints, is this really a
runtime parameter rather than a topology parameter? If so, the integrator
will have a reference to the constraint scheme and could configure it.

Implicitly created objects (e.g. the neighborlist and electrostatics)
have parameters that can be documented with the rest of the class,
so we should not use the matplotlib strategy of passing ``**kwargs`` along
such that it is hard to figure out what options are available and how they
are processed. If we want to provide convenience, we could bundle options
as dictionaries to be passed to, e.g., an nlist_params argument to the
integrator.

.. sourcecode :: python

     # With no arguments, use current gmx code to detect and allocate compute resources.
     context = gmx.Context()
     # With default 'context=None', implicitly use gmx.Context()
     minimization = gmx.System(context=context, structure=atoms, topology=topology)
     minimizer = gmx.md.Steep(emtol=500., nsteps=1000, coulombtype='PME', nlist='grid')

     # In practice, kwargs will likely come from parameters files.
     # classes or module attibutes may have shorthand string names.
     minimization_params = {'emtol': 500,
                            'nsteps': 1000}
     integration_params = {'coulombtype': 'PME',
                           'constraints': [],
                           'nlist': 'grid'}
     nlist_params = {'frequency': 1,
                     'rcut': 1.0}

     # bind the Integrator.
     minimization.integrator(minimizer)

For implicit creation and binding, use the class
methods ``gmx.System.from_*()``, which avoid overly-complicated System
constructor and can be extended to package common sets of parameters
and simple workflows. E.g. maybe ``gmx.System.from_minimization()``
Maybe the ``from_`` is cumbersome.

Optionally, add loggers

.. sourcecode :: python

     # optionally
     # get an energy group of all atoms in the system
     egroup = gmx.Group(minimization.all_atoms())
     minimization.reporters.append(gmx.reporter.LogEnergy(period=1, groups=[egroup]))


Run energy minimizer.
e.g. ``mdrun -v -deffnm em-vac``

.. sourcecode :: python

     # optionally set output behavior
     minimization.filename_base('em-vac')
     # Run with execution context implicitly configured
     minimization.run()

Solvate.
~~~~~~~~

Fill the box with water. e.g.
``genbox -cp em-vac.gro -cs spc216.gro -p fws.top -o fws-b4ion.gro``

After adding atoms, we will use the same integration method with the
same electrostatics, topology, and neighborlist parameters, but a few
simulation parameters change.

For user-friendliness, AtomData can use getattr(x, to_pdata) or
something to see if an automatic conversion is possible, or objects in the
gmx.solvent submodule could already be AtomData objects.
Similarly, the solute and solvent arguments in solvate() could try to cast
to AtomData objects.

.. sourcecode :: python

     solvent = gmx.AtomData('spc216.gro') # load solvent molecule coordinates
     grofile = gmx.AtomData('em-vac.gro') # load energy-minimized configuration
     atoms = gmx.util.solvate(solute=grofile,
                      solvent=solvent,
                      topology=gmx.Topology('fws.top'))

     # or maybe
     solvent = gmx.AtomData(gmx.data.spc216)
     # and we are still holding the topology object, which was extended
     # earlier with a solvent definition from gmx.amber99sb-ildn.tip3p_flexible
     atoms = gmx.utils.solvate(solute=grofile, solvent=solvent, topology=topology)

     # or (alternate use)

     atoms = minimization.atoms # load energy-minimized configuration
     atoms = gmx.util.solvate(configuration=atoms,
                      solvent=solvent,
                      topology=topology)
     # Note that 'atoms' now refers to a new object and will need to be reattached.
     system.load_configuration(atoms)

     # or (implementation option)

     # Use utility to solvate a loaded system
     atoms = gmx.util.solvate(system=minimization)

     # Optionally save configuration to file
     atoms.save('fws-b4ion.gro')

Note that 'atoms' now refers to a new object and will need to be reattached with
``system.load_configuration(atoms)``
This invalidates the neighborlist, domain decomposition, etc.
The validity of the topology (i.e. does it define the solvent?) could be
checked with the solvate command, but should definitely be checked whenever
load_configuration() is called. If done at the call, we impose a requirement
that topology must be updated before configuration, which seems reasonable.
Again, some refinement may still need to occur conceptually on the
encapsulation of structure data versus topology data in terms of configuration,
atom typing, molecule / residue type definitions, full system topology,
and miscellaneous metadata, such as the "names" and such used for file
I/O and/or due to differences in conventions for the contents of molecular
data files.

Add ions to solvated system.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

e.g.
``grompp -f em-sol-pme.mdp -c fws-b4ion.gro -p fws.top -o ion.tpr``

Presumably we need a tpr to get the PME and neighbor parameters from the mdp file?
Need to figure out what is really needed by genion. It may be that it is
more appropriate as an integrator, like Steep(). Maybe we need a different
term. Integrator is too specific, and so is MD. Simply Updater?

.. sourcecode :: python

     minimization = gmx.System.from_mdpfile(
                                     mdpfile='em-sol-pmd.mdp',
                                     topology=topology,
                                     atoms=gmx.AtomData('fws-b4ion.gro'),
                                     )
     minimization.save('ion.tpr')

Use utility to insert ions. e.g. ``genion -s ion.tpr -o fws-b4em.gro -neutral -conc 0.15 -p fws.top -g ion.log``

Solvate, insert-molecule, and genion perform similar functions with different
algorithms and options.
This should probably be
reconsidered as something more abstract. I.e. an alchemy module or
add/change atom methods to invoke these algorithms with appropriate
parameters as arguments.

.. sourcecode :: python

     # Should the utility be allowed to modify the input in place?
     # e.g. gmx.util.genion(system=minimization)
     # For early iterations, I think not...
     # Create a new System from the input System
     minimization = gmx.util.genion(system=minimization, conc=0.15, neutral=True)

     # implementation alternative: return AtomData

     atoms = gmx.util.genion(system=minimization, conc=0.15, neutral=True)
     minimization.load_configuration(atoms)

     # or (implementation alternative) is genion sufficiently coupled
     # to System objects to simply be a fancy construction helper?

     minimization = gmx.System.genion(system=minimization, conc=0.15, neutral=True)

     # Alternatively...

     minimization.change_residues(mode='solvate', ...)
     minimization.change_residues(mode='pack', ...) # genion
     # mode='trans' or mode='trans_rot'
     # for the two filling mechanims in insert-molecule

     minimization.save('fws-b4em.gro')

Minimize energy in solvated system.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

e.g. ``grompp -f em-sol-pme.mdp -c fws-b4em.gro -p fws.top -o em-sol.tpr``

.. sourcecode :: python

     # Prep simulation
     minimization = gmx.System.from_mdpfile(
                                     mdpfile='em-sol-pmd.mdp',
                                     topology=topology,
                                     atoms=gmx.AtomData('fws-b4em.gro'),
                                     )

     # Optionally: suppress output configured in mdp file
     for r in minimization.reporters:
         minimization.reporters.remove(r) # or maybe `del r`

     # Optionally, If we already have a handle to mimimizer, we can just reuse it.
     minimimizer.set_param(emtol=250.0, nsteps=5000)
     # The minimizer will notice that its convergence is invalidated both by the
     # new value of emtol and by the updated atoms. The new nsteps parameter
     # is used the next time it tries to converge.

     # mdrun -v -deffnm em-sol
     minimization.filename_base('em-sol')
     minimization.run() # generates trajectory and other configured output

Two step equilibration.
~~~~~~~~~~~~~~~~~~~~~~~

Compare to
::

     $ grompp -f nvt-pr-md.mdp -c em-sol.gro -p fws.top -o nvt-pr.tpr
     $ mdrun -deffnm nvt-pr
     $ grompp -f npt-pr-md.mdp -c em-sol.gro -p fws.top -o npt-pr.tpr
     $ mdrun -deffnm npt-pr

Set up and run two System objects in sequence.

.. sourcecode :: python

     # grompp -f nvt-pr-md.mdp -c em-sol.gro -p fws.top -o nvt-pr.tpr
     # mdrun -deffnm nvt-pr
     nvt = gmx.System.from_mdpfile(mdpfile='nvt-pr-md.mdp',
                                  topology=topology,
                                  atoms='em-sol.gro',
                                  fname_base='nvt-pr')
     # implementation option: allow other parameter value types
     nvt = gmx.System.from_mdpfile(mdpfile='nvt-pr-md.mdp',
                                  topology=topology,
                                  atoms=minimization,
                                  fname_base='nvt-pr')
     # optionally save input record
     nvt.save('nvt-pr.tpr')

     nvt.run()

     # grompp -f npt-pr-md.mdp -c em-sol.gro -p fws.top -o npt-pr.tpr
     # mdrun -deffnm npt-pr
     npt = gmx.System.from_mdpfile(mdpfile='npt-pr-md.mdp',
                                  topology=topology,
                                  atoms='em-sol.gro',
                                  fname_base='npt-pr')
     # or
     npt = gmx.System.from_mdpfile(mdpfile='npt-pr-md.mdp',
                                  topology=topology,
                                  atoms=nvt,
                                  fname_base='npt-pr')

     npt.save('npt-pr.tpr')
     npt.run()

Construct topology and structure objects
----------------------------------------
There remains a question of the relationships between the AtomData,
Topology (system topology), and Molecule (molecular topology) structures.
I think this will be sorted out in part through examination of current Gromacs
architecture and in part through user input and/or inspiration from popular
tools like MDTraj and MDAnalysis. The interface to these data structures may
or may not have much in common with the implementation, as there is a lot of
opportunity for relational data and lazy attributes.

Note that system topologies and structure data are likely to be generated at the
same time, but common use cases may involve substituting entire sets of coordinates
or system topologies (in part or in whole).

From pdb and force field definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

From raw data and force field definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

From a currently instantiated simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

From generic structure and force field information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

From bare metal?
~~~~~~~~~~~~~~~~

TODO

Run with frozen N and C terminals
---------------------------------

Make index for residue groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

Configure simulation with frozen terminals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e.g.

| ``make_ndx –f clg_b4md.pdb –o clg_ter.ndx``
| Prompted input: ``r 1-36 & a C N CA`` for residue selection, get the new group number and rename for convenience: ``name 15 Terminal``, then ``v`` to view and verify.
| Additional MDP entries:
|    ``energygrp_excl = Terminal Terminal Terminal SOL``
|    ``freezegrps = Terminal``
|    ``freezedim = Y Y Y``
| ``grompp –f md.mdp –c pr.gro –p clg.top –n clg_ter.ndx –o md.tpr``

Run simulation with index
~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

Trajectory Manipulation tasks
-----------------------------

Extract a simulation frame
``trjconv -f traj.xtc -s file.tpr -o time_3000ps.pdb -dump 3000``

Access atom coordinates from the trajectory

| 1a. Create Trajectory object from file
| 1b. Get handle to Trajectory or Frame object held by simulation
| 1c. Use Trajectory object created via API
| 2. Retrieve Frame

Re-center a molecule
``trjconv –f traj.xtc –o traj_center.xtc –s str_b4md.gro –pbc nojump -center``

Make the dummy gro file for the g_covar analysis.
``trjconv –s ../md.tpr –f dangle.trr –o resiz.gro –n covar.ndx –e 0``

Concatenate trajectories
``trjcat –f md1.xtc md2.xtc md3.xtc ... (etc) –o mdall.xtc -settime``

Trajectory analysis
-------------------

The following example is from a comment by Teemu on Redmine issue `1625 <https://redmine.gromacs.org/issues/1625#note-13>`_ with changes reflecting the following considerations.

1. add instances of analyzer classes instead of keywords
2. call-back at level of Runner can be supplemented if custom modules can easily be subclassed from AnalysisModule
3. Data flow management is a a project aim and will be available Pythonically.

::

     # Adapted from https://redmine.gromacs.org/issues/1625#note-13
     import gmx.analysis

     class MyCustomDataCombiner(gmx.analysis.Module):
         def __init__(self, distances, angles):
             gmx.analysis.Module.__init__() # does the base class do anything?
             self._distances = distances
             self._angles = angles

         def process_frame(self):
             curr_dist = self._distances.get_current_frame()
             curr_angles = self._angles.get_current_frame()
             # do whatever custom analysis on the combined angles and distances

     # It is useful to require that operations be atomically added to a runner
     # after the execution context has been initialized so that each node
     # can be initialized once (at binding) and have the opportunity to raise exceptions.
     # Create execution context, detect environment, and parse launch parameters
     runner = gmx.Runner(**kwargs)

     # Parameters can specify the trajectory, begin and end times etc
     # Lazy initialization defers any I/O until object is run in the Gromacs execution context.
     traj = gmx.Trajectory.from_file(filename=name, **kwargs)
     # traj uses a special Trajectory class method to create an object with no externally-accessible input connection.

     # Create data flow graph by adding the first node
     runner.add_module(traj)
     # traj is initialized and should do as much error checking as possible.

     # These make the runner run the specified predetermined analysis modules, with the given parameters used to initialize them.
     # Bind two inputs that Distance will consume when run.
     # Some type checking or other sanity tests can be performed here.
     distances = gmx.analysis.Distance(traj, traj)
     # Further configure the modular task by setting parameters.
     # Missing parameters will be hard to detect until the graph is run or
     # a runner (implicit or explicit) initializes the module, presumably when it is mapped into the execution context.
     distances.params(select=['<selection1>', '<selection2>'], **kwargs)

     # bind Analyzer object to execution context
     runner.add_module(distances)
     # distances is initialized and can raise exceptions for, e.g., missing parameters, insufficient/inappropriate resources.

     angles = runner.add_module(gmx.analysis.Angle(g1='vector', g2='vector', group1=[...], group2=[...]))

     # The custom module gets called each frame, and can combine the data from the various
     # intermediate data structures that the predefined modules provide to compute, e.g.,
     # cross-correlation.
     custom = MyCustomDataCombiner(distances, angles)

     runner.add_frame_callback(custom.process_frame)

     # This reads in the trajectory
     # and does all the specified analyses with a single pass.
     runner.run()

In the above example, graph nodes are bound to their upstream data source at construction.

Alternatively, binding inputs could be a separate step, allowing more flexibility for nodes at the data flow terminus. E.g.::

     # Create Distances object
     distances = gmx.analysis.Distances()
     distances.params(select=['<selection1>', '<selection2>'], **kwargs)

     # Create empty Trajectory object
     traj = gmx.Trajectory()

     # Over-ride dangling input connection with special method
     traj.from_file(filename, **kwargs)

     # Bind distances input
     distances.input(traj, traj)

     runner.add_module(traj)
     runner.add_module(distance)
     # Note the ultimate output of the graph may be another Trajectory object...


It might be important for additional methods to handle fancier processing or multiple outputs for some modules.
::

     trajectory = Trajectory(filename, mode="r")

     # Distances is an Analyzer subclass that requires two input objects that provide an iterator for Nx3 arrays or an iterator of Frame objects, or something of the sort.
     distances1 = Distances()
     distances1.input(trajectory(selection=group1), trajectory(selection=group2)) # connect data sources
     distances2 = Distances()
     distances2.input(trajectory(selection=group1), trajectory(selection=group3))

     # CrossCorrelation is an Analyzer subclass that requires two input objects that provide an iterator for Nx3 arrays or of Vector objects or something else, depending on design choices.
     xcorr = CrossCorrelation()
     xcorr.input(distances1.output(), distances2.output()) # connect data sources

The binding operation input() could throw a ValueError or TypeError if the inputs do not provide objects of the right type.
For Analyzers that provide multiple outputs, the subscriber could either ask for a specific standardized attribute (like 'traj') or require a specific binding, either through attributes or more granular objects.  E.g.::

     bar.input(foo.traj)

Analysis tools
=====================

For more flexible workflows, whether driven from python or some other (TBD) external API,
the line will be blurred between analysis tools, utilities such as ``trajconv``,
and actual simulation.

Early proof-of-concept targets for abstraction and API access will be
I/O, file parsing, manipulations such as ``trajconv`` components, and the more actively developed analysis modules using the new C++ trajectory analysis API.
We expect to be influenced by dataflow-centric prior art such as the TensorFlow framework.

Whether explicit or implicit, the gmx framework will include some sort of "runner" that manages the graph of operations,
exploiting data localization in the Gromacs parallelization schemes on behalf of the user.



Read-only access internal data structures
=========================================
::

    md_params = gmx.read_mdp(filename)
    context = gmx.initialize(**md_params)
    system = gmx.Integrator(context)
    system.run(1e6)
    print(system.force_rec.get_nonbonded().get_potential_energy())
    system.run(1e6)
    system.dump()
    context.finalize()

Alternatively or additionally, a slot for a Reporter or Analyzer class could be
inserted into the simulation loop.

OpenMM-like::

    system.reporters.append(ThermoReporter('data.csv', kineticEnergy=True, potentialEnergy=True, pressure=True))
    # low performance Python-level custom reporters?
    class MyReporter(Reporter):
        def report():
            cur_fr = self.system._cpp_integrator.force_rec
            sys.stderr.write(cur_fr.get_nonbonded().get_potential_energy())
    system.reporters.append(MyReporter())

HOOMD-like::

    context.initialize()
    ...
    system = init_from_xml()
    ...
    # hook into a loggable quantity provided by the Lennard-Jones force-calculation kernel
    logger = analyze.log('out.dat', period=1, quantities=['potential_energy_lj'])
    system.run(10) # 'out.dat' has ten lines
    e = logger.query('potential_energy_lj') # retrieve most recently logged value

Retrievable handles to modules owned by MD engine
=================================================
::

    coordinates, topology = gmx.utils.pdb2gmx(**{**ff_params, **water_params, ...})
    minimization = gmx.integrator.Steep(coords=coordinates, topology=topology, **params)
    ff = minimization.get_force_field()
    minimization.run(nsteps)
    for compute in ff.get_computes():
        print(compute.name, compute.get_energy())

Modules exposed at a "context" level
====================================
such that a handle exists before and after the MD integrator runs.

.. sourcecode :: python

    with gmx.initialize_empty(**runtime_params) as context:
        coordinates, topology = gmx.utils.pdb2gmx(
                                        force_field=gmx.forcefield.amber99sb-ildn,
                                        pdb_file="fws.pdb",
                                        coords_file="fws.gro",
                                        topology_file="fws.top",
                                        water=gmx.watermodel.tip3p)
        context.set_positions(coordinates)
        integrator = VVIntegrator(topology, ff, **md_params)
        context.set_integrator(integrator)
        context.reporters.append(gmx.reporter.EnergyLogger(**log_params))
        # integrator, force compute kernels, and reporter are not yet initialized

        context.run(1)
        # domain decomposition has taken place, particles migrated, compute

        # modules initialized for the topology, and all internal data structures initialized
        context.step(1000)
        # nothing changed due to user input, so nothing was reinitialized

        context.nonbonded[0] # handle to the first force kernel applied during nonbonded evaluation
        context.nonbonded.append(Plumed.FancyForce(...))
        context.step(1)
        # if cut-off range has changed, domain decomposition has been updated, etc.

        # Now get energy from last compute, i.e. Plumed.FancyForce
        energy = context.nonbonded[-1].get_energy()
        # communication has been triggered and all ranks now have the total energy stored locally in Python.

        # Interact with the distributed objects at a low level.
        positions = gmx.get_comm().Gather(context.get_local_particles().positions[:])

        print("Done")

    # 'with' calls context.__exit__() to make sure context.finalize() cleans up
    # despite remaining references

TensorFlow analogy
==================
The following is a mind-expanding thought experiment, not a personal goal.

.. sourcecode :: python

    coordinates, topology = gmx.utils.pdb2gmx(force_field=gmx.forcefield.amber99sb-ildn,
                                    pdb_file="fws.pdb",
                                    coords_file="fws.gro",
                                    topology_file="fws.top",
                                    water=gmx.watermodel.tip3p)
    # Create timestep variable
    current_step = gmx.Step(0)

    # Create particle positions variable
    xyz = gmx.ParticleCoordinates()
    xyz.initialize(coordinates) # operation deferred until "graph" is evaluated

    # Create velocities variable
    v = gmx.Velocities()
    v.initialize(gmx.ThermalVelocities(**params))

    # Create a node to hold system state at a time step
    state = gmx.SystemState()

    # Create operation to invoke non-bonded force compute kernel
    nonbonded = gmx.ForceCompute(xyz, topology)

    # Create operation to invoke constraints kernel
    constraints = gmx.Constraints(xyz, topology)

    # Create operation to invoke electrostatics kernel
    pme = gmx.Electrostatics(xyz, **params)

    # Operation to sum forces for integrator
    forces = gmx.SumForces(nonbonded, constraints, pme)

    # Operation to perform integration, updating data handles and state
    state = gmx.Integrate(xyz, v, forces)

    # Operation to perform simulation loop to advance to a given time step
    output = gmx.RunUpTo(state, current_step, gmx.Step(1000))

    with gmx.initialize(**runtime_params) as session:
        # Create an operation that updates to the next timestep
        for _ in range(10):
            session.eval(output)
            print("Step: {}, electrostatic energy: {}".format(
                    session.eval(current_step),
                    session.eval(pme).energy))
