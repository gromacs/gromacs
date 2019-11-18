The modular simulator {#page_modularsimulator}
==============================================

A new modular approach to the GROMACS simulator is described. The
simulator in GROMACS is the object which carries out a simulation. The
simulator object is created and owned by the runner object, which is
outside of the scope of this new approach, and will hence not be further
described. The simulator object provides access to some generally used
data, most of which is owned by the runner object.

## Using the modular simulator
GROMACS will automatically use the modular simulator for the velocity
verlet integrator (`integrator = md-vv`), if the functionality chosen
in the other input parameters is implemented in the new framework.
Currently, this includes NVE simulations, NVT simulations (
`tcoupl = v-rescale` only), NPH simulation (`pcoupl = parrinello-rahman` 
only), and NPT simulations (`tcoupl = v-rescale` and 
`pcoupl = parrinello-rahman` only), with or without free energy perturbation.

To disable the modular simulator for cases defaulting to the new framework,
the environment variable `GMX_DISABLE_MODULAR_SIMULATOR=ON` can be set. To
use the new framework also for `integrator = md` (where the functionality is
implemented), the environment variable `GMX_USE_MODULAR_SIMULATOR=ON` can 
be set to override legacy default.

## Legacy implementation

In the legacy implementation, the simulator consisted of a number of
independent functions carrying out different type of simulations, such
as `do_md` (MD simulations), `do_cg` and `do_steep` (minimization),
`do_rerun` (force and energy evaluation of simulation trajectories),
`do_mimic` (MiMiC QM/MM simulations), `do_nm` (normal mode analysis),
and `do_tpi` (test-particle insertion).

The legacy approach has some obvious drawbacks:
* *Data management:* Each of the `do_*` functions defines local data,
  including complex objects encapsulating some data and functionality,
  but also data structures effectively used as "global variables" for
  communication between different parts of the simulation. Neither the
  ownership nor the access rights (except for `const` qualifiers) are
  clearly defined.
* *Dependencies:* Many function calls in the `do_*` functions are
  dependent on each others, i.e. rely on being called in a specific
  order, but these dependencies are not clearly defined.
* *Branches:* The flow of the `do_*` functions are hard to understand
  due to branching. At setup time, and then at every step of the
  simulation run, a number of booleans are set (e.g. `bNS` (do neighbor
  searching), `bCalcEner` (calculate energies), `do_ene` (write
  energies), `bEner` (energy calculation needed), etc). These booleans
  enable or disable branches of the code (for the current step or the
  entire run), mostly encoded as `if(...)` statements in the main `do_*`
  loop, but also in functions called from there.
* *Task scheduling:* Poorly defined dependencies and per-step branching
  make task scheduling (e.g. parallel execution of independent tasks)
  very difficult.
* *Error-prone for developers:* Poorly defined dependencies and unclear
  code flow make changing the simulator functions very error-prone,
  rendering the implementation of new methods tedious.

## The modular simulator approach

The main design goals of the new, fully modular simulator approach
include
* *Extensibility:* We want to ease maintenance and the implementation
  of new integrator schemes.
* *Monte Carlo:* We want to add MC capability, which can be mixed with
  MD to create hybrid MC/MD schemes.
* *Data locality & interfaces:* We aim at localizing data in objects,
  and offer interfaces if access from other objects is needed.
* *Multi-stepping:* We aim at a design which intrinsically supports
  multi-step integrators, e.g. having force calls at different
  frequencies, or avoid having branches including rare events
  (trajectory writing, neighbor search, ...) in the computation loops.
* *Task parallelism:* Although not first priority, we want to have a
  design which can be extended to allow for task parallelism.

The general design approach is that of a **task scheduler**. *Tasks*
are argument-less functions which perform a part of the computation.
Periodically during the simulation, the scheduler builds a
*queue of tasks*, i.e. a list of tasks which is then run through in
order. Over time, with data dependencies clearly defined, this
approach can be modified to have independent tasks run in parallel.

The approach is most easily displayed using some pseudo code:

    class ModularSimulator : public ISimulator
    {
        public:
            //! Run the simulator
            void run() override;
        private:
            std::vector<ISignaller*> signallers_;
            std::vector<ISimulatorElement*> elements_;
            std::queue<SimulatorRunFunction*> taskQueue_;
    }

    void ModularSimulator::run()
    {
        constructElementsAndSignallers();
        setupAllElements();
        while (not lastStep)
        {
            // Fill the task queue with new tasks (can be precomputed for many steps)
            populateTaskQueue();
            // Now simply loop through the queue and run one task after the next
            for (auto task : taskQueue)
            {
                (*task)();  // run task
            }
        }
    }

This allows for an important division of tasks.

* `constructElementsAndSignallers()` is responsible to **store the
  elements in the right order**. This includes the different order of
  element in different algorithms (e.g. leap-frog vs. velocity
  verlet), but also logical dependencies (energy output after compute
  globals).
* `populateTaskQueue()` is responsible to **decide if elements need to
  run at a specific time step**. The elements get called in order, and
  decide whether they need to run at a specific step. This can be
  pre-computed for multiple steps. In the current implementation, the
  tasks are pre-computed for the entire life-time of the neighbor
  list.
* **Running the actual simulation tasks** is done after the task queue
  was filled.  This is achieved by simply looping over the task list,
  no conditionals or branching needed.

### Simulator elements

The task scheduler holds a list of *simulator elements*, defined by
the `ISimulatorElement` interface. These elements have a
`scheduleTask(Step, Time)` function, which gets called by the task
scheduler. This allows the simulator element to register one (or more)
function pointers to be run at that specific `(Step, Time)`. From the
point of view of the element, it is important to note that the
computation will not be carried out immediately, but that it will be
called later during the actual (partial) simulation run. From the
point of view of the builder of the task scheduler, it is important to
note that the order of the elements determines the order in which
computation is performed. The task scheduler periodically loops over
its list of elements, builds a queue of function pointers to run, and
returns this list of tasks. As an example, a possible application
would be to build a new queue after each domain-decomposition (DD) /
neighbor-searching (NS) step, which might occur every 100 steps. The
scheduler would loop repeatedly over all its elements, with elements
like the trajectory-writing element registering for only one or no
step at all, the energy-calculation element registering for every
tenth step, and the force, position / velocity propagation, and
constraining algorithms registering for every step. The result would
be a (long) queue of function pointers including all computations
needed until the next DD / NS step, which can be run without any
branching.

### Signallers

Some elements might require computations by other elements. If for
example, the trajectory writing is an element independent from the
energy-calculation element, it needs to signal to the energy element
that it is about to write a trajectory, and that the energy element
should be ready for that (i.e. perform an energy calculation in the
upcoming step). This requirement, which replaces the boolean branching
in the current implementation, is fulfilled by a Signaller - Client
model. Classes implementing the `ISignaller` interface get called
*before* every loop of the element list, and can inform registered
clients about things happening during that step. The trajectory
element, for example, can tell the energy element that it will write
to trajectory at the end of this step. The energy element can then
register an energy calculation during that step, being ready to write
to trajectory when requested.

### Sequence diagrams

#### Pre-loop
In the loop preparation, the signallers and elements are created and
stored in the right order. The signallers and elements can then
perform any setup operations needed.

\msc
hscale="2";

ModularSimulator,
Signallers [label="ModularSimulator::\nSignallers"],
Elements [label="ModularSimulator::\nElements"],
TaskQueue [label="ModularSimulator::\nTaskQueue"];

--- [ label = "constructElementsAndSignallers()" ];
    ModularSimulator => Signallers [ label = "Create signallers\nand order them" ];
    ModularSimulator => Elements [ label = "Create elements\nand order them" ];
--- [ label = "constructElementsAndSignallers()" ];
|||;
|||;

--- [ label = "setupAllElements()" ];
    ModularSimulator => Signallers [ label = "Call setup()" ];
    Signallers box Signallers [ label = "for signaler in Signallers\n    signaller->setup()" ];
    |||;
    ModularSimulator => Elements [ label = "Call setup()" ];
    Elements box Elements [ label = "for element in Elements\n    element->setup()" ];
--- [ label = "setupAllElements()" ];
\endmsc

#### Main loop
The main loop consists of two parts which are alternately run until the
simulation stop criterion is met. The first part is the population of
the task queue, which determines all tasks that will have to run to
simulate the system for a given time period. In the current implementation,
the scheduling period is set equal to the lifetime of the neighborlist.
Once the tasks have been predetermined, the simulator runs them in order.
This is the actual simulation computation, which can now run without any
branching.

\msc
hscale="2";

ModularSimulator,
Signallers [label="ModularSimulator::\nSignallers"],
Elements [label="ModularSimulator::\nElements"],
TaskQueue [label="ModularSimulator::\nTaskQueue"];

ModularSimulator box TaskQueue [ label = "loop: while(not lastStep)" ];
ModularSimulator note TaskQueue [ label = "The task queue is empty. The simulation state is at step N.", textbgcolor="yellow" ];
|||;
|||;
ModularSimulator box ModularSimulator [ label = "populateTaskQueue()" ];
ModularSimulator =>> TaskQueue [ label = "Fill task queue with tasks until next neighbor-searching step" ];
|||;
|||;
ModularSimulator note TaskQueue [ label = "The task queue now holds all tasks needed to move the simulation from step N to step N + nstlist. The simulation for these steps has not been performed yet, however. The simulation state is hence still at step N.", textbgcolor="yellow" ];
|||;
|||;

ModularSimulator => TaskQueue [ label = "Run all tasks in TaskQueue" ];
TaskQueue box TaskQueue [label = "for task in TaskQueue\n    run task" ];
TaskQueue note TaskQueue [ label = "All simulation computations are happening in this loop!", textbgcolor="yellow" ];
|||;
|||;
ModularSimulator note TaskQueue [ label = "The task queue is now empty. The simulation state is at step N + nstlist.", textbgcolor="yellow" ];
ModularSimulator box TaskQueue [ label = "end loop: while(not lastStep)" ];

\endmsc

#### Task scheduling
A part of the main loop, the task scheduling in `populateTaskQueue()` 
allows the elements to push tasks to the task queue. For every scheduling 
step, the signallers are run first to give the elements information about 
the upcoming scheduling step. The scheduling routine elements are then 
called in order, allowing the elements to register their respective tasks.

\msc
hscale="2";

ModularSimulator,
Signallers [label="ModularSimulator::\nSignallers"],
Elements [label="ModularSimulator::\nElements"],
TaskQueue [label="ModularSimulator::\nTaskQueue"];

--- [ label = "populateTaskQueue()" ];
    ModularSimulator box ModularSimulator [ label = "doDomainDecomposition()\ndoPmeLoadBalancing()" ];
    ModularSimulator =>> Elements [ label = "Update state and topology" ];
    |||;
    |||;

    ModularSimulator note ModularSimulator [ label = "schedulingStep == N\nsimulationStep == N", textbgcolor="yellow" ];
    ModularSimulator box TaskQueue [ label = "loop: while(not nextNeighborSearchingStep)" ];
        ModularSimulator => Signallers [ label = "Run signallers for schedulingStep" ];
        Signallers box Signallers [label = "for signaller in Signallers\n    signaller->signal(scheduleStep)" ];
        Signallers =>> Elements [ label = "notify" ];
        Signallers note Elements [ label = "The elements now know if schedulingStep has anything special happening, e.g. neighbor searching, log writing, trajectory writing, ...", textbgcolor="yellow" ];
        |||;
        |||;

        ModularSimulator => Elements [ label = "Schedule run functions for schedulingStep" ];
        Elements box Elements [label = "for element in Elements\n    element->scheduleTask(scheduleStep)" ];
        Elements =>> TaskQueue [ label = "Push task" ];
        Elements note TaskQueue [ label = "The elements have now registered everything they will need to do for schedulingStep.", textbgcolor="yellow" ];
        ModularSimulator note ModularSimulator [ label = "schedulingStep++", textbgcolor="yellow" ];

    ModularSimulator box TaskQueue [ label = "end loop: while(not nextNeighborSearchingStep)" ];
--- [ label = "populateTaskQueue()" ];
ModularSimulator note ModularSimulator [ label = "schedulingStep == N + nstlist\nsimulationStep == N", textbgcolor="yellow" ];

\endmsc

## Acceptance tests and further plans

In January 2019, we defined acceptance tests which need to be 
fulfilled to make the modular simulator the default code path:
* End-to-end tests pass on both `do_md` and the new loop in
  Jenkins pre- and post-submit matrices
* Physical validation cases pass on the new loop
* Performance on different sized benchmark cases, x86 CPU-only
  and NVIDIA GPU are at most 1% slower -
  https://github.com/ptmerz/gmxbenchmark has been developed to
  this purpose.

After the MD bare minimum, we will want to add support for
* Pulling
* Full support of GPU (current implementation does not support
GPU update)

Using the new modular simulator framework, we will then explore
adding new functionality to GROMACS, including
* Monte Carlo barostat
* hybrid MC/MD schemes
* multiple-time-stepping integration

We will also explore optimization opportunities, including
* re-use of the same queue if conditions created by user input are 
  sufficiently favorable (by design or when observed)
* simultaneous execution of independent tasks

We will probably not prioritize support for (and might consider
deprecating from do_md for GROMACS 2020)
* Simulated annealing
* REMD
* Simulated tempering
* Multi-sim
* Membrane embedding
* QM/MM
* FEP lambda vectors
* Fancy mdp options for FEP output
* MTTK
* Essential dynamics
* Constant acceleration groups
* Ensemble-averaged restraints
* Time-averaged restraints
* Freeze, deform, cos-acceleration

## Signallers and elements

The current implementation of the modular simulator consists of
the following signallers and elements:

### Signallers

All signallers have a list of pointers to clients, objects that
implement a respective interface and get notified of events the
signaller is communicating.

* `NeighborSearchSignaller`: Informs its clients whether the
  current step is a neighbor-searching step.
* `LastStepSignaller`: Informs its clients when the current step
  is the last step of the simulation.
* `LoggingSignaller`: Informs its clients whether output to the
  log file is written in the current step.
* `EnergySignaller`: Informs its clients about energy related
  special steps, namely energy calculation steps, virial
  calculation steps, and free energy calculation steps.
* `TrajectoryElement`: Informs its clients if writing to
  trajectory (state [x/v/f] and/or energy) is planned for the
  current step. Note that the `TrajectoryElement` is not a
  pure signaller, but also implements the `ISimulatorElement`
  interface (see section "Simulator Elements" below).

### Simulator Elements

#### `TrajectoryElement`
The `TrajectoryElement` is a special element, as it
is both implementing the `ISimulatorElement` and the `ISignaller`
interfaces. During the signaller phase, it is signalling its
_signaller clients_ that the trajectory will be written at the
end of the current step. During the simulator run phase, it is
calling its _trajectory clients_ (which do not necessarily need
to be identical with the signaller clients), passing them a valid
output pointer and letting them write to trajectory. Unlike the
legacy implementation, the trajectory element itself knows nothing
about the data that is written to file - it is only responsible
to inform clients about trajectory steps, and providing a valid
file pointer to the objects that need to write to trajectory.

#### `StatePropagatorData`
The `StatePropagatorData` takes part in the simulator run, as it might
have to save a valid state at the right moment during the
integration. Placing the StatePropagatorData correctly is for now the
duty of the simulator builder - this might be automated later
if we have enough meta-data of the variables (i.e., if
`StatePropagatorData` knows at which time the variables currently are,
and can decide when a valid state (full-time step of all
variables) is reached. The `StatePropagatorData` is also a client of
both the trajectory signaller and writer - it will save a
state for later writeout during the simulator step if it
knows that trajectory writing will occur later in the step,
and it knows how to write to file given a file pointer by
the `TrajectoryElement`.

#### `EnergyElement`
The `EnergyElement` takes part in the simulator run, as it
does either add data (at energy calculation steps), or
record a non-calculation step (all other steps). It is the
responsibility of the simulator builder to ensure that the
`EnergyElement` is called at a point of the simulator run
at which it has access to a valid energy state.

#### `ComputeGlobalsElement`
The `ComputeGlobalsElement` encapsulates the legacy calls to
`compute_globals`. While a new approach to the global reduction
operations has been discussed, it is currently not part of this
effort. This element therefore aims at offering an interface
to the legacy implementation which is compatible with the new
simulator approach.

The element currently comes in 3 (templated) flavors: the leap-frog case,
the first call during a velocity-verlet integrator, and the second call
during a velocity-verlet integrator. It is the responsibility of the
simulator builder to place them at the right place of the
integration algorithm.

#### `ForceElement` and `ShellFCElement`
The `ForceElement` and the `ShellFCElement` encapsulate the legacy
calls to `do_force` and `do_shellfc`, respectively. It is the
responsibility of the simulator builder to place them at the right
place of the integration algorithm. Moving forward, a version of these
elements which would allow to only calculate forces of subsets of
degrees of freedom would be desirable to pave the way towards multiple
time step integrators, allowing to integrate slower degrees of freedom
at a different frequency than faster degrees of freedom.

#### `ConstraintElement`
The constraint element is implemented for the two cases of constraining
both positions and velocities, and only velocities. It does not change the constraint
implementation itself, but replaces the legacy `constrain_coordinates`
and `constrain_velocities` calls from update.h by elements implementing
the ISimulatorElement interface and using the new data management.

#### `Propagator`
The propagator element can, through templating, cover the different
propagation types used in the currently implemented MD schemes. The
combination of templating, static functions, and having only the
inner-most operations in the static functions allows to have performance
comparable to fused update elements while keeping easily re-orderable
single instructions.

Currently, the (templated) implementation covers four cases:
 * *PositionsOnly:* Moves the position vector by the given time step
 * *VelocitiesOnly:* Moves the velocity vector by the given time step
 * *LeapFrog:* A manual fusion of the previous two propagators
 * *VelocityVerletPositionsAndVelocities:* A manual fusion of VelocitiesOnly 
    and PositionsOnly, where VelocitiesOnly is only propagated by half the 
    time step of PositionsOnly.

The propagators also allow to implement temperature and pressure coupling
schemes by offering (templated) scaling of the velocities.

#### `CompositeSimulatorElement`
The composite simulator element takes a list of elements and implements
the ISimulatorElement interface, making a group of elements effectively
behave as one. This simplifies building algorithms.

#### `VRescaleThermostat`
The `VRescaleThermostat` implements the v-rescale thermostat. It takes a
callback to the propagator and updates the velocity scaling factor
according to the v-rescale thermostat formalism.

#### `ParrinelloRahmanBarostat`
The `ParrinelloRahmanBarostat` implements the Parrinello-Rahman barostat.
It integrates the Parrinello-Rahman box velocity equations, takes a
callback to the propagator to update the velocity scaling factor, and
scales the box and the positions of the system.

#### `FreeEnergyPerturbationElement`
The `FreeEnergyPerturbationElement` holds the lambda vector and the
current FEP state, offering access to its values via getter
functions. The FreeEnergyPerturbationElement does update the lambda
values during the simulation run if lambda is non-static. It
implements the checkpointing client interface to save its current
state for restart.

## Data structures

### `StatePropagatorData`
The `StatePropagatorData` contains a little more than the pure
statistical-physical micro state, namely the positions,
velocities, forces, and box matrix, as well as a backup of
the positions and box of the last time step. While it takes
part in the simulator loop to be able to backup positions /
boxes and save the current state if needed, it's main purpose
is to offer access to its data via getter methods. All elements
reading or writing to this data need a pointer to the
`StatePropagatorData` and need to request their data explicitly. This
will later simplify the understanding of data dependencies
between elements.

Note that the `StatePropagatorData` can be converted to and from the
legacy `t_state` object. This is useful when dealing with
functionality which has not yet been adapted to use the new
data approach - of the elements currently implemented, only
domain decomposition, PME load balancing, and the initial
constraining are using this.

### `EnergyElement`
The EnergyElement owns the EnergyObject, and is hence responsible
for saving energy data and writing it to trajectory. It also owns
the tensors for the different virials and the pressure as well as
the total dipole vector.

It subscribes to the trajectory signaller, the energy signaller,
and the logging signaller to know when an energy calculation is
needed and when a non-recording step is enough. The simulator
builder is responsible to place the element in a location at
which a valid energy state is available. The EnergyElement is
also a subscriber to the trajectory writer element, as it is
responsible to write energy data to trajectory.

The EnergyElement offers an interface to add virial contributions,
but also allows access to the raw pointers to tensor data, the
dipole vector, and the legacy energy data structures.

### `TopologyHolder`
The topology object owns the local topology and holds a constant reference
to the global topology owned by the ISimulator.

The local topology is only infrequently changed if domain decomposition is
on, and never otherwise. The topology holder therefore offers elements to register
as ITopologyHolderClients. If they do so, they get a handle to the updated local 
topology whenever it is changed, and can rely that their handle is valid 
until the next update. The domain decomposition element is defined as friend 
class to be able to update the local topology when needed.

The topology holder is not a `ISimulatorElement`, i.e. it does not take part in the
simulator loop.

## Infrastructure
### `DomDecHelper` and `PmeLoadBalanceHelper`
These infrastructure elements are responsible for domain decomposition and 
PME load balancing, respectively. They encapsulate function calls which are 
important for performance, but outside the scope of this effort. They rely 
on legacy data structures for the state (both) and the topology (domdec).
    
The elements do not implement the ISimulatorElement interface, as
the Simulator is calling them explicitly between task queue population
steps. This allows elements to receive the new topology / state before
deciding what functionality they need to run.

### `Checkpointing`
The `CheckpointHelper` is responsible to write checkpoints. In the
longer term, it will also be responsible to read checkpoints, but this
is not yet implemented.

Writing checkpoints is done just before neighbor-searching (NS) steps,
or before the last step. Checkpointing occurs periodically (by default,
every 15 minutes), and needs two NS steps to take effect - on the first
NS step, the checkpoint helper on master rank signals to all other ranks
that checkpointing is about to occur. At the next NS step, the checkpoint
is written. On the last step, checkpointing happens immediately before the
step (no signalling). To be able to react to last step being signalled,
the CheckpointHelper does also implement the `ISimulatorElement` interface,
but does only register a function if the last step has been called.

Checkpointing happens at the top of a simulation step, which gives a
straightforward re-entry point at the top of the simulator loop.

In the current implementation, the clients of CheckpointHelper fill a
legacy t_state object (passed via pointer) with whatever data they need
to store. The CheckpointHelper then writes the t_state object to file.
This is an intermediate state of the code, as the long-term plan is for
modules to read and write from a checkpoint file directly, without the
need for a central object. The current implementation allows, however,
to define clearly which modules take part in checkpointing, while using
the current infrastructure for reading and writing to checkpoint.
