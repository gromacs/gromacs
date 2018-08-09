=======================
External C++ Interfaces
=======================

A library and headers are under development to allow stable C++ API access to
|Gromacs| installations. Progress is tracked at Redmine issue
`2585 <https://redmine.gromacs.org/issues/2585>`_ and sub-issues.

gmxapi Milestones
=================

Project milestones are identified with "gmxapi milestone #", where "#" is an
integer that will not change during the project, once declared. Milestones
roughly correspond to current or future Redmine issues in scope.
Large scale design dependencies are illustrated in the accompanying chart.
Each numbered feature (chart node) is expected to be from 1 to 10 Gerrit
changes, generally 3 to 5.

.. graphviz::

    digraph milestones{
        graph [rankdir = "LR"];

        gmx1  [label="gmx1  Design documentation strategy /\n project management plan", shape=box, style=rounded];

        gmx2  [label="gmx2  Versioned libgmxapi target \n for build, install, headers, docs", shape=box,style=rounded];
        gmx1 -> gmx2

        gmx3  [label="gmx3  Integration testing", shape=box,style=rounded];
        gmx2 -> gmx3

        gmx4  [label="gmx4  Library access to \nMD runner symbols", shape=box,style=rounded];
        gmx2 -> gmx4

        gmx5  [label="gmx5  Provide runner with context manager\n from which to get factory functions", shape=box,style=rounded];
        gmx4 -> gmx5

        gmx6  [label="gmx6  Extensible MDModules and ForceProviders", shape=box,style=rounded];
        gmx1 -> gmx6

        gmx7  [label="gmx7  Binding API for higher-level client code", shape=box,style=rounded];
        gmx6 ->gmx7
        gmx5 -> gmx7

        gmx8  [label="gmx8  Binding API for plug-in ForceProviders", shape=box,style=rounded];
        gmx2 -> gmx8

        gmx9  [label="gmx9  Headers and adapter classes\n for Restraint framework", shape=box,style=rounded];
        gmx8 -> gmx9

        gmx10 [label="gmx10 MD signalling API", shape=box,style=rounded];
        gmx5 -> gmx10
        gmx23 -> gmx10

        gmx11 [label="gmx11 Replace MPI_COMM_WORLD\n with a Context resource", shape=box,style=rounded];
        gmx5 -> gmx11
        gmx12 -> gmx11

        gmx12 [label="gmx12 Runtime API for \n sharing / discovering hardware / parallelism resources", shape=box,style=rounded];
        gmx2 -> gmx12

        gmx13 [label="gmx13 API for working directory, input/output targets?", shape=box,style=rounded];
        gmx12 -> gmx13

        gmx14 [label="gmx14 Generalized pull groups / generalized sites", shape=box,style=rounded];
        gmx4 -> gmx14

        gmx15 [label="gmx15 API logging resource", shape=box,style=rounded];
        gmx5 -> gmx15

        gmx16 [label="gmx16 Exception handler firewall", shape=box,style=rounded];
        gmx4 -> gmx16
        gmx17 -> gmx16

        gmx17 [label="gmx17 API status object", shape=box,style=rounded];
        gmx1 -> gmx17

        gmx18 [label="gmx18 Thinner test harness", shape=box,style=rounded];
        gmx19 -> gmx18
        gmx20 -> gmx18
        gmx21 -> gmx18

        gmx19 [label="gmx19 API manipulation\n of simulation input / output", shape=box,style=rounded];
        gmx1 -> gmx19

        gmx20 [label="gmx20 Accessible test data resources", shape=box,style=rounded];

        gmx21 [label="gmx21 Break up runner state\n into a hierarchy of RIAA classes\n with API hooks", shape=box,style=rounded];
        gmx5 -> gmx21

        gmx22 [label="gmx22 API management of input objects
         gmx22.1 Structure, topology
         gmx22.2 Microstate
         gmx22.3 Simulation state
         gmx22.4 Simulation parameters
         gmx22.5 Runtime parameters / execution environment
         gmx22.x Anything else?", shape=box,style=rounded];
         gmx21 -> gmx22
         gmx19 -> gmx22

        gmx23 [label="gmx23 Event hooks or signals
          gmx23.1 checkpoint
          gmx23.2 time step number
          gmx23.3 input configuration
          gmx23.4 input topology
          gmx23.5 input state
          gmx23.6 simulation parameters
          gmx23.7 output data streams", shape=box,style=rounded];
          gmx4 -> gmx23

        gmx24 [label="gmx24 API expression of\n MDOptions interfaces\n and embedded user documentation", shape=box,style=rounded];
        gmx7 -> gmx24

        gmx25 [label="gmx25 Replace std::exit (gmx_fatal)\n with exceptions", shape=box, style=rounded];
        gmx1 -> gmx25

        gmx26 [label="gmx26 API messaging resources", shape=box, style=rounded];
        gmx5 -> gmx26

    #    gmx27 [label="gmx27 Integration testing", shape=box, style=rounded];
    #    gmx2 -> gmx27

        gmx28 [label="gmx28 set simulation parameters from API", shape=box, style=rounded];
        gmx4 -> gmx28

        gmx29 [label="gmx29 grompp functionality", shape=box, style=rounded];
        gmx1 -> gmx29

        gmx30 [label="gmx30 manipulate files / topologies", shape=box, style=rounded];
        gmx1 -> gmx30
    }


gmxapi milestone 1: Design documentation strategy / project management plan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Refer to https://redmine.gromacs.org/issues/2585

gmxapi milestone 2 (Issue #2586) Versioned libgmxapi target for build, install, headers, docs

Establish a framework in which further API infrastructure can be developed.
Installations of GROMACS 2019 and on should provide a build environment that is
forward compatible with client software developed against newer features to allow
straight-forward incremental feature development and integration testing.

Reference Redmine issue `2586 <https://redmine.gromacs.org/issues/2586>`_.

gmxapi milestone 3 Integration testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Gmxapi interfaces should continue functioning with unchanged semantics for other GROMACS changes, or API level needs to be incremented according to semantic versioning.
* External projects need to be tested outside of the gromacs build tree to validate external interfaces of installation. Suggested external projects: Python package, sample_restraint, yet-to-be-written integration test suite.
* Tests should be clear about the API version they are testing, and we should test all versions that aren’t unsupported (though we need a policy in this regard) and we can note whether new API updates are backwards compatible.
* Forward-compatibility testing: we should at least _know_ whether we are breaking old client code and include release notes, regardless of policy
* ABI compatibility testing? (should we test mismatched compilers and such?)
* Example code in documentation should be tested, if possible.

gmxapi milestone 4 Library access to MD runner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* mdrun CLI program is an API client
* Non-CLI client code built against gmxapi can launch MD simulations.

Relates to #2229
Relates to https://redmine.gromacs.org/issues/2375

gmxapi milestone 5 Provide runner with context manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reference Redmine issue `2587 <https://redmine.gromacs.org/issues/2587>`_.

gmxapi milestone 6 Extensible MDModules and ForceProviders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ForceProviders obtained after tMPI threads have spawned.
* MDModules list extended at runtime during simulation launch.
* External code may be provided to the runner to instantiate or get a handle to a module.
* Expanded Context class can broker object binding by registering and holding factory functions for modules, as well as other resources that may be implemented differently in different environments.
* Somewhere in here, MDModules either need access to the integral timestep number or the ability to register call-backs or signals on a schedule.

Relates to #2590, #2574, #1972

Do MDModules live in a scope of tight association with an integrator? Do we need other concepts, like RunnerModules? Or subdivisions like MDForceModule, MDObserverModule, MDControlModule?

gmxapi milestone 7 Binding API for higher-level client code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gmxapi milestone 8 Binding API for plug-in ForceProviders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ultimately tied to gmx5 and gmx24, but we can start stabilizing the external interfaces now. The external interfaces
are for (a) user interface / workflow management code, and (b) MD extension code. We define a simple message-passing
C structure along with PyCapsule name and semantics. An MD extension object can provide a factory method with which
the MD Runner can get an IMDModules interface at simulation launch. The object pointed to may exist before and/or
after the lifetime of the simulation. It must be understood that the IMDModule handle will be obtained on every rank.
Design should consider future infrastructure and needs, but does not need to implement now. (expressing data
dependencies and locality, negotiating parallelism, expressing periodicity) Short-term implementation may require
workarounds for some of these, but the workaround can mostly be segregated from this issue’s resolution.

Relates to #2590

gmxapi milestone 9 Headers and adapter classes for Restraint framework
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Relates to #1972, #2590, https://github.com/kassonlab/sample_restraint

gmxapi milestone 10 MD signalling API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Relates to #2224

gmxapi milestone 11 Replace MPI_COMM_WORLD with a Context resource
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Part of reducing dependence on global variables and definitions. Allow client
code to define the MPI group in which a simulation operates.

gmxapi milestone 12 Runtime API for hardware and parallelism resources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Runtime API for sharing / discovering hardware / parallelism resources

* Libgmxapi requests resources from libgromacs from the current node
* CUDA environment can be manipulated but we shouldn’t have to deal with that for a while
* Evolving task scheduling interfaces, expressing data locality
* Concepts of time and timestep

gmxapi milestone 13 API for working directory, input/output targets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unify / clarify the API for the various "context"-related resources.

gmxapi milestone 14 Generalized pull groups / “generalized sites”
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Continue to develop facilities like LocalAtomSet with facilities to negotiate
local data availability.

gmxapi milestone 15 API logging resource
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Log "file" artifacts are produced through API, allowing extensibility and abstraction from filesystem dependence. Progress has already been made in this direction, but the logging resource could be more clearly owned by the client code (or a Context object owned or managed on behalf of the client code) rather than created and destroyed in, say, the Mdrunner.

gmxapi milestone 16 Exception handler firewall
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

currently the gmx binary has a commandline runner thing that catches the exceptions, reports an error and exits, but the API can and should do something else, because it plays the same role as the commandline runner

gmxapi milestone 17 API status object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Status type defines the interface for discovering operation success or failure, plus details.
* Consistent status object interface is portable across Python, C++, and C
* Status object can be used to ferry information across API boundaries from exceptions thrown. Exceptions could be chained / status nested.

Questions:

* What are concerns and solutions for memory allocation for status objects? Should objects own one or generate one on function return?
* Should the API (or Context) keep a Status singleton? A Status stack? Or should operations create ephemeral Status objects, or objects implementing a Status interface?
* Should the status object contain strings, reference strings mapped by enum, or defer textual messages to messaging and logging facilities?

gmxapi milestone 18 Thinner test harness
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Much of the code related to testing sets up harnesses for command line tools.
This is often not helpful for testing API-driven simulations, but should also
be replaced by access to standardized API facilities.

gmxapi milestone 19 API manipulation of simulation input / output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Also aids testing.

GlobalTopology class and IGlobalTopologyUser interface underway will help here, so that
client changes to the global topology can ripple through to the modules because the ones that care have registered themselves at setup time

gmxapi milestone 20 Accessible test data resources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data for testing mdrun is not accessible to other test suites.

gmxapi milestone 21 Break up runner with a clear set of states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mdrunner is reentrant in that it launches itself for tMPI. There are not clear
points at which user input has been handled, parallel resources have been set up,
and the simulation launches. Break up into a hierarchy of RAII classes with
API hooks.

* break up mdrun program into clearly defined layers and phases
* CLI program parses various inputs in order to launch an Mdrunner object that is CLI-agnostic
* launching tMPI threads and other significant changes of state establish a sequence or hierarchy of invariants through RAII and/or State pattern.
* Sebastian Wingbermuhle working now on aspects of this for hybrid MC/MD (ref #2375, …)

Closely related to gmxapi milestone 4

gmxapi milestone 22 API management of input objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Structure, topology
* Microstate
* Simulation state
* Simulation parameters
* Runtime parameters / execution environment
* Anything else?

gmxapi milestone 23 Event hooks or signals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Event hooks or signals for

* checkpoint
* time step number or delta / trajectory advancement
* input configuration
* input topology
* input state
* simulation parameters
* output data streams

gmxapi milestone 24 API expression of MDOptions interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Module developers are able to express MDP options, command line options,
and embedded user documentation in a fairly compartmentalized way in for CLI
mdrun. We need to figure out how that maps to a run-time extensible API and a
Builder pattern in which user interface (like command-line argument processing)
is handled in Director code before implementation objects have been created.

gmxapi milestone 25 Avoid sys::exit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generally, replace std::exit (gmx_fatal)with exceptions

* Root out gmx_fatal, clearly define regular exit points and exception throwers
* API firewall should catch exceptions from gmx and convert to status objects for ABI compatibility. (gmx17)
* Clearly document regular and irregular shutdown behavior under MPI, tMPI, and generally, specifying responsibilities
* Create issue tickets for discovered missing exception safety, memory leaks, opportunities for RAII refactoring, and complicated protocols that should either be better documented or replaced with a clearer hierarchy (or sequence) of invariants

gmxapi milestone 26 API messaging resources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Abstraction for status messages, such as are currently printed to stdout or stderr

gmxapi milestone 27 (retracted)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gmxapi milestone 28 set simulation parameters from API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Short term: mdrun CLI-like functionality to override other input is sufficient

Long term: sufficient API to update parameters between phases of simulation work

Implementation roadmap is probably

1. Inject argv fields
2. Write to input_rec or other structures
3. Interact with MDOptions framework

gmxapi milestone 29 API access to grompp functionality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Generate runnable input from user input
* United implementation for workflow API and utility functions (e.g. possibility of deferred execution / data transfer)
* Ultimately should not require writing output to (tpr) file
* File inputs ultimately should be generalized to API objects

gmxapi milestone 30 API access to file and data manipulation tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* United implementation for workflow API and utility functions (e.g. possibility of deferred execution / data transfer)
* Utility API should be sufficient to reimplement CLI tools
* I/O should ultimately be separate from algorithm; filesystem interaction optional
* Consider feature requirements of other projects such as MDAnalysis.

Scope
=====

There are definitely design points for consideration that are left out of this list merely because they are not essential to gmxapi functionality or because gmxapi doesn’t have strong dependence on the ultimate design choice. These topics include:

* Task scheduling framework
* Insertion points in the MD loop
* Encapsulation of integrator

Further downstream, this infrastructure is necessary to support new high level interfaces to GROMACS, but the discussion of such interfaces is deferred as much as possible to separate issues to streamline incorporation of the changes proposed here in less public / stable code.
