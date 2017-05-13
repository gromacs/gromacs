=================
Execution context
=================

From low-level parallelism to high-level task scheduling and data flow,
implementation details for features vary in ways that can be hidden from the
high-level interfaces, depending on the compute resources.
The execution context is compartmentalized with an abstraction layer.

* Hide the details of work scheduling to allow separation of compute resource
  allocation from task preparation
* Provide a layer in which work scheduling and data flow can be managed
  by external software or extended with minimally invasive code additions.
* Allow object-oriented interaction with handles to objects that are not well defined outside
  of a running context (*e.g.* object-oriented interfaces to internal modules,
  a node in a data graph that is not yet executing, a plugin in a simulation
  that has either finished or not yet started)

Communicator objects should be represented and have appropriate interfaces at
each API layer that is aware of the execution context.

Extending the facilities of the context abstraction is a long term project, but
should be considered in short-term design. Most immediately, it will guide the
handling of code paths for single-threaded versus multithreaded and MPI features.

For instance, it should be easy to write a Python script that can run a Gromacs
workflow on a single core or in an MPI context, but it should also be easy for
a user to make sensible use of their available concurrency APIs. An API client
(Python or C++) should be allowed, but not required, to adapt code execution
for data locality.

It should be difficult and explicit for API clients to
invalidate distributed data structures or inadvertently trigger expensive
communication, but appropriate notifiers should be available so that such
operations do not violate encapsulation (e.g. signal/slot, publisher-subscriber,
observer). The same framework may be used to announce and handle data events in
graph execution.
