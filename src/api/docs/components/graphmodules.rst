==================
Runner and modules
==================

We envision a convergence of the module managers and runners from
the simulation and analysis frameworks. Ultimately, we expect a
data flow manager to be able to connect and manage modules generating and
operating on particle configurations and trajectories including
modules that implement simulation algorithms to generate trajectories.

Part of this vision includes at least one additional abstraction layer
between the module manager and the user interface, unlike the current
architecture in which the modules and runners are deeply coupled to
their command line invocation.

Modular tasks for simulation and analysis can be implemented in libgromacs or external software to be executed in a Gromacs execution context with a common runner (see :doc:`context`). Development iterations will improve functionality and performance of pipelined simulation/analysis and more complex data graphs. High-level user-friendly tools can be implemented in terms of API features.

Trivial access to low-level modular functionality (e.g. evaluating energy of a configuration) may be allowed through a simpler API component (see :doc:`microstate`) and/or a distinct trivial module runner.
