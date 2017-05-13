=============================
Architecture and design notes
=============================

..  important concepts, design principles,
    architecture-level requirements, protocols and interfaces
    (sequence diagrams, object diagrams, class diagrams (where
    available))

Principles, requirements, and guidelines.

Principles
==========

Performance versus Flexibility
------------------------------

Make it easy for the user to do efficient operations with data and constructs owned by the framework; make it easy but transparent and explicit when the user performs flexible but less efficient operations moving data between Gromacs ownership and external ownership.

API data ownership
------------------
Crossing the boundary of data ownership by the framework versus the user interface should be explicit.

As a consequence of the above, data "owned" by the graph of gmx modules can be optimized for performance and minimal communication. Crossing that boundary by getting data to and from the python interpreter should be explicit. Implicit casting of data to GMX objects should not be performed if it crosses this boundary. Fetching data should require an explicit and syntactically distinct accessor method.

Optimize tool chaining within the Gromacs framework
---------------------------------------------------

Similarly (and motivated by some of the explicit data movement patterns above), we will want to expose some primitives for variables under Gromacs management (simple array operations) as well as simplify and modularize some common operations performed in analysis workflows between Gromacs tasks.  If the common usage is GMX.operationA -> X -> GMX.operationB and X is common or sufficiently simple, then it should be possible to do X in the GMX framework.  If X does not occur between two GMX operations, X is highly complex, and/or X is not commonly used, then the case for management by the GMX framework as opposed to external code is less strong.

Catch errors early
------------------
Gromacs tools are often structured the way they are so that the user has the
chance to get errors before submitting a job to a queue or using resources.
We should keep this in mind. To every extent possible, we should make sure that
the runner will not fail due to configuration errors when it is ready to schedule
work.

Protocols and Interfaces
------------------------
The API must be sufficiently
self-documenting to facilitate correct use *and* extension. Required behavior
of clients and derived classes, and expected behavior of API objects, must be
clear for stateful interfaces, *e.g.* protocols requiring a sequence of multiple
API calls or API calls with side effects. To help ensure that documented design
matches actual behavior, design documentation should be embedded in and extracted
from the implementation, tagged with API version metadata, and verified with regression tests.

This is not a port or a wrapper
-------------------------------
Don't expose old cruft and unnecessary functionality.

One example that comes up repeatedly: aspects of pdb2gmx and trajconv only exist for historic reasons and need not be exposed, duplicated, or preserved. Other aspects should be discrete API calls. We should maintain an eye for tasks better provided by other tools and aim for interoperability. We should provide an easy and well-documented mechanism for bringing outside tools into the Gromacs environment where appropriate for performance. I.e. there may be a gmx.analysis.map() mechanism or some such to apply arbitrary computation under the management of Gromacs workflow and (parallel) resource management.

High level tools and convenience APIs are to be built on a C++ API accessible
from outside of the core library. This public API is built with the library at a
lower level than any UI-aware code. This requires reworking aspects of the library
API in ways that are also helpful to the maintenance of current code and tools.

While CLI and API use cases are distinct, CLI tools should be implemented at the
highest possible layer to ensure uniform semantics and maintainability for common
functionality. Many current tools have CLI or terminal-centric implementation at
a low level. Adding public API functionality is an opportunity to refactor such
tools, but some fundamental idioms will need to be revisited, such as interactive
options refinement.

Required or included software
=============================

Suggestions:

* Python >= 3.4
* pybind11
* Eigen

Validation and verifiability
============================

The Gromacs project rightly emphasizes the reproducibility of scientific results
and the ability to guarantee that simulation results can be tied to a known,
validated software set. These priorities are reflected in design decisions and
cultivated paradigms related to the simulation input and output files.

This documentation discusses what records and
bookkeeping exist in API-driven Gromacs use to accomodate this scientific due
diligence.

Specifically, what is the role of the content of a script, itself? What is
sufficient to describe for posterity the parameters of a specific set of output?
What discussion, documentation, clarification, etc, if any, needs to happen
regarding mutable input parameters, as in ensemble simulations or adaptive
workflows? Is there a canonical way to record parameters necessary for
reproducibility or simply a set of conventions or standards? If there is not a
canonical representation (e.g. TPR file), how and to what degree are tools
implemented for automated sanity checking, etc.?

Requirements
------------

Simulation output should be clearly mapped to the inputs and code that produced
it.

Model, protocol, and method should be easily and concisely captured for
reproducibility of published (publishable) work.

There should be no additional ambiguity introduced by API use. E.g. The meaning
and implications of a TPR file should not change, or should only change as a
result of the broader Gromacs roadmap.
