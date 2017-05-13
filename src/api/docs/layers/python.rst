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
