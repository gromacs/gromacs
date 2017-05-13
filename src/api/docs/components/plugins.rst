====================
Simulation internals
====================

Simulation plugins, call-backs, and API-facilitated introspection or manipulation
require some additional facilities and revised data structures.

* Hooks in the simulation loop for commonly extended behaviors to allow non-invasive extensibility through plugins
* Publicly derivable interface classes for extension code not (yet) included with Gromacs
* run-time binding of plugins via API
* granular access to internal machinery for call-backs, rapid prototyping, etc.

Some considerations follow.

Objects / Data owned by the integrator are limited to data
generated and used specifically because of the integration
method or algorithm and used to update the configuration
from one simulation step to the next.

* thermostat/barostat internal state
* constraints?
* virial?

Objects used by the Integrator

* simulation box definition
* PRNG
* force / potential evaluators
* particle and topology data
* logger
* messenger
* communicator(s)
* load balancer(s)
* user-provided parameters
* neighbor list
* particle selection / grouping structures
* checkpointing facilities

The integrator, and objects and functors used by the Integrator, should be created
at a higher scope than the simulation loop so that they can
provide facilities to other code before and after simulation
as well as during simulation to objects bound before the loop
begins. Public facilities can include facilities provided to the
integrator and/or functionality and data of interest to users
or other objects.

* applied forces, last calculated values
* contribution to energy and virial
* current parameter values
* state introspection
* parallel / distributed access methods as well as gather operations
* useful "signals" or publish/subscribe hooks for data owned

Some basic abstract interfaces will be obvious, while more
esoteric code should be free to provide extended interfaces.

Less public functionality may still be exposed, possibly through
additional wrappers or a sequence of calls to account for inner loop
optimizations. E.g. ``update()`` or ``calculate()`` methods called by
the simulation runner every step, or methods to invalidate some
aspect of simulation state to force recalculation or reinitialization
of something.

Methods to invalidate state should only be called by the object that
owns the relevant data, since whether or not a certain action invalidates
state may depend on implementation. E.g. changing the size of the simulation
box doesn't invalidate domain decomposition if there is only one domain,
while rescaling particle positions invalidates Verlet/neighbor lists for
truncated non-bonded potentials, but not necessarily all force compute
kernels.

In the absence of the Boost signals code, I have seen the
[nano-signal-slot](https://github.com/NoAvailableAlias/nano-signal-slot) package used effectively
to provide signals and slots for handling things like, *e.g.*
`BoxChanged`, `ParticleNumChanged`...
Initial implementation will require only a few signals to connect core modules
(for *e.g.* triggering domain decomposition or neighbor data rebuilds)
but that additional signals and handlers will be added as needed,
and that the facility will be very useful for
developers of extension code, such as alchemical techniques (`ParticleTypesChanged`,
`ChargesUpdated`, `NonIntegratorParticleTranslation`, `NBRangeChanged`) and
encapsulated workflows (`TSetChanged`, `NDegreesFreedomChanged`).
