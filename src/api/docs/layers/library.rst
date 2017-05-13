====================
Core Gromacs library
====================

Plugins
=======

It should be possible for users to implement external potentials, particle
interactions, live simulation analysis, and high-performance analysis-driven
simulation tweaking without altering the distributed Gromacs code by writing
to a plugin API. Target use cases are tools like Plumed, experimental free energy
analysis, long-term code development that may or may not eventually be submitted
for inclusion, and everyone who has
hacked on the external force part of mdrun to implement custom pulling or biasing
methods.

In addition to the obvious basic access to internal data structures, it will be
helpful to expand the functionality of the communicator(s) to allow local code
to request data to be sent or retrieved so that plugins can take advantage of
the existing parallelism without incurring additional communication. To support
future task-based parallelism, additional infrastructure may be necessary for
managed data access handles and mechanisms for selecting/grouping atoms.

Plugins should also be provided with facilities for checkpointing, logging,
status messages, and whatever else they need to avoid doing their own I/O.

I suspect there is even less burden of ABI compatibility for plugins, but the
public API, library API, and plugin API should be have separately managed
versioning to facilitate API compatibility checks.

The API should provide plugins with (at least) the following resources:

* logging facility
* message facility
* access to communicators
* ability to request local data for next step
* iterators for particles in whatever structural data types exist
* serialization / checkpoint hooks
* neighbor list(s)

Callback hooks
==============

It should be as easy as possible to create C++ objects that can
be plugged into the simulation loop as viewers, observers, actors,
or just compute kernels, but a simple callback API allows greater
flexibility. But how should a callback be bound and invoked?

I propose one call-back function can be registered with the runner to be called
every `period` steps (provide a helper function for time units).
The function will be called with a single argument that is a
reference to the simulation object.
The call-back function can access whatever statefull information
it needs through the API and/or be configured with necessary
references to other simulation objects before the loop is run.

Supporting infrastructure
=========================

Updates in the core Gromacs library are necessary to support compartmentalization, but are distinct from this project. Such dependencies belong in an issue tracking system, but can be linked from here for convenience before public release.

Various abstraction in libgromacs will be necessary in these projects. Some is already underway.

* Terminal I/O
* File I/O abstraction

    * trajectory output
    * input and configuration files
    * log files

* Prompted user interaction
* CLI and runtime options
* Simulation configuration, input state
* Gathered parallel data, running state
* task-level modules
* Compute kernels used within ``do_md()`` (for plug-ins as well as external access)
* Domain--rank binding / Interprocess communication / MPI communicator
