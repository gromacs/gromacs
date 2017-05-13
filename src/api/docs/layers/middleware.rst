=============================
Execution context abstraction
=============================

Concurrency API / Data flow API / compute resource management.

Many useful tools already provide Python APIs for integration with compute
resources, data flow managers, and workflow managers. The abstractions planned
for the Python interface should lend themselves to optimal interaction when
invoked from such resources, but we can also think about adding a layer between
the user interface and the execution context to allow workflows to be executed
with minimal alteration on a workstation, HPC cluster, pilot manager, or cloud
computing resource.

We may be able to take advantage of existing code for data flow management /
data graph execution, but the high level APIs should be able to rely on deferring
execution of chains (or graphs) of API calls to allow for optimizations to be
made based on available resources and data locality.

Gromacs should have its own workflow checkpointing system.
