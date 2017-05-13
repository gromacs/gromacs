======
Guides
======

Build and install
=================

How to
======

..
    links to documentation by example begin as user stories.
    Cross-link with components doc (feature) when possible.
    Unavailable workflows belong in a scrum board or issue tracking system, but
    are too noisy for the main Gromacs Redmine. For the moment, curate them here
    to clarify targeted features.

As a Python user
----------------

Use bottled workflows / tasks
so that I can get going quickly.

Access coordinates from a trajectory file
so that I can use my favorite Python-based analysis tool.

Pass non-file-backed trajectory data to and from the MD engine
so that I can avoid filesystem overhead.

Transform a PDB record to a model suitable for a Gromacs input record
so that I can use molecular data from a common source.

Run -> analyze -> run in a pipeline
so that filesystem and MPI overhead is avoided

Sanity-check my input parameters
so that I can avoid precision or other numerical errors or known bad physics leading to invalid simulation results (e.g. as provided by grompp)


As an advanced workflow
-----------------------
Have flexible ways to represent/pass data BLOBs and streams
so that I can adapt to native APIs for containers, workflow managers, and other cyberinfrastructure.

Deal with proxy objects
so that I don't trigger unnecessary communications across successive API calls.

Coordinate batches of simulation results or trajectory streams
so that I can asynchronously configure and launch the next simulation ASAP.

Define data flow or graphs for deferred execution
so that I can let lower level software optimize execution and communications.

Rely on robust checkpointing of data graph state
so that I can avoid duplicated computation after network interruptions or job failure.

Include simulation and analysis tasks in the same data graph
so that I can use a single, simpler interface.

As an API client
----------------
Manipulate simulation box and atom configuration
so that I can solvate, adjust ions, or set periodic boundary conditions

Separate core module parameters from user interface so that I can process options
without additional API calls.

Make a set of well-defined changes to a template of simulation input
so that I can systematically prepare the next batch of (uncoupled) simulations.

Compute the potential energy for a provided atom configuration
so that I can analyze a single configuration easily.

Measure the potential energy of a given Selection and force field during a simulation
so that I can compare force fields or perform analyses

Retrieve microstate data from a live (or paused) simulation
so that I can assure the mapping between calculations performed within and outside of Gromacs regardless of trajectory file format conventions, restrictions, or overhead.

Produce/edit a topology through a high-level interface
so that I can use optimized Gromacs kernels without constraint to a compiled force field implementation.

Use a concurrency API in conjunction with the Gromacs API
so that I can take advantage of parallelism and data locality optimizations from versatile external code.

As a plugin
-----------

Apply an external field to an atom selection
so that I can implement additional physics.

Define new potentials by implementing a simple interface
so that I can access native Gromacs parallelism and generalizable optimizations.

Access objects passed in by the user
so that workflow-scope parameters can be managed in a single place.

Access contributions of individual force field components (potentials) to calculated values
so that model details are accessible with consistent semantics.

Request, access, and announce data through an interface to the library communicator(s)
so that I can have the local data I need without the ability to trigger communication.
