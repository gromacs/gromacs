mdrun modules {#page_mdmodules}
=============

All functionality of mdrun that is not basic computation of bonded and
non-bonded forces and integration should be implemented through the
mdrun "MD module" interface. This is a limited interface that completely
hides the internals of the mdrun machinery. The intent is to keep this
interface simple and stable. Thus most of the "special" functionalities,
both those are part of the \Gromacs main code base and developed
by users, keep functioning with each new version of \Gromacs.

Many current MD modules are "force providers" implemented in a single call-back
function that is called at each MD step during the force calculation.
Here coordinates and parameters of local atoms are provided and forces,
the energy and virial can be returned. Additionally there are several
callbacks during setup and domain decomposition time, and interfaces for
declaring MDP file options and writing output. This enables all
the specifics of the setup and input parameter handling to be contained
inside the MD module. Extensions to support other kinds of operations are
envisaged and proposals welcome.

To keep the API limited, the functionality is limited, also performance
wise. There is e.g. no formal support for computation on GPUs.
Note that any MPI communication of an MD module can significantly
impact the parallel performance of mdrun, as mdrun avoids and bundles
global communication for standard runs as much as possible.
However, when most of the standard force computation is handled by GPUs,
the MD modules computation and communication can overlap with the standard
work by using the, otherwise empty, CPU resources.

Several current MD modules are in sub-directories in the
`src/gromacs/applied-forces/` directory (note that not all modules
there use the MD module interface yet). There is a, very simple,
electric field module, as well as a more complex neural network
potential module, interfaces to the Colvars and Plumed packages, and
the density-fitting module.
When implementing a new module, it can be useful to have a look
at these modules.

The rest of this page documents those parts of the modularity mechanism that
have taken a clear form.  Generalizing and designing other parts may require
more code to be converted to modules to have clearer requirements on what the
mechanism needs to support and what is the best way to express that in a
generally usable form.

Structure of a module
---------------------

Each module implements a factory that returns an instance of gmx::IMDModule.
This interface has methods that in turn refer to other interfaces:
gmx::IMdpOptionProvider, gmx::IMDOutputProvider, and gmx::IForceProvider.
The module also implements these interfaces (or a subset of them), and code
outside the module only calls methods in these interfaces.

See documentation of the individual interfaces for details of what they
support.

Implementation of a module
--------------------------

Modules are constructed by composition of interfaces (i.e. abstract classes,
general with pure virtual methods lacking implementations), so that e.g.
trajectory-writing code can loop over containers of pointers to
gmx::IMDOutputProvider without needing to know about all the concrete types
that might implement that interface.

The module classes should not be extended by using them as a base
class, which is expressed with the final keyword in the class
definition. Generally, modules will implement different flavours of
functionality, perhaps based on user choices, or available computing
resources. This should generally be implemented by providing variable
behavior for the methods that are called through the above
interfaces. Either code should branch at run time upon some data
contained by the module (e.g. read from the mdp options), or that the
module class should contain a pointer to an internal interface class
whose concrete type might be chosen during setup from the set of
implementations of that internal interface. Such an approach keeps
separate the set of interfaces characteristic of "MD modules" from
those that are particular to flavours of any specific module.

The virtual methods that the module classes inherit from their
interfaces should be declared as `override`, to express the intent
that they implement a virtual function from the interface. This
permits the compiler to check that this is true, e.g. if the interface
class changes. The `virtual` keyword should not be specified,
because this is redundant when `override` is used. This follows
the Cpp Core Guidelines (guideline C.128).

Handling mdp input
------------------

To accept parameters from an mdp file, a module needs to implement
gmx::IMdpOptionProvider.

initMdpOptions() should declare the required input parameters using the options
module.  In most cases, the parameters should be declared as nested sections
instead of a flat set of options.  The structure used should be such that in
the future, we can get the input values from a structured mdp file (e.g., JSON
or XML), where the structure matches the declared options.  As with other uses
of the options module, the module needs to declare local variables where the
values from the options will be assigned.  The defined structure will also be
used for storing in the tpr file (see below).

initMdpTransform() should declare the mapping from current flat mdp format to
the structured format defined in initMdpOptions().  For now, this makes it
possible to have an internal forward-looking structured representation while
the input is still a flat list of values, but in the future it also allows
supporting both formats side-by-side as long as that is deemed necessary.

On the implementation side, the framework (and other code that interacts with
the modules) will do the following things to make mdp input work:

* When grompp reads the mdp file, it will first construct a flat
  KeyValueTreeObject, where each input option is set as a property.

  It then calls initMdpTransform() for the module(s), and uses the produced
  transform to convert the flat tree into a structured tree, performing any
  defined conversions in the process.  This transformation is one-way only,
  although the framework keeps track of the origin of each value to provide
  sensible error messages that have the original mdp option name included.

  It calls initMdpOptions() for the module(s), initializing a single Options
  object that has the input options.

  It processes the structured tree using the options in two steps:

  * For any option that is not specified in the input, it adds a property to
    the tree with a default value.  For options specified in the input, the
    values in the tree are converted to native values for the options (e.g.,
    from string to int for integer options).
  * It assigns the values from the tree to the Options object.  This will make
    the values available in the local variables the module defined in
    initMdpOptions().

  Note that currently, the module(s) cannot use storeIsSet() in options to know
  whether a particular option has been provided from the mdp file.  This will
  always return true for all the options.  This is a limitation in the current
  implementation, but it also, in part, enforces that the mdp file written out
  by `gmx grompp -po` cannot produce different behavior because of set/not-set
  differences.

* grompp -po writes an mdp file that was equivalent to the input,
  which is implemented by calling buildMdpOutput() for each module, to
  prepare a builder object that is used with writeKeyValueTreeAsMdp().
  As with the old flat tree, the values given by the user's input are
  preserved, but not the ordering of options, or their formatting.

* When grompp writes the tpr file, it writes the structured tree (after the
  default value and native value conversion) into the tpr file.

* When mdrun reads the tpr file, it reads the structured tree.
  It then broadcasts the structure to all ranks.  Each rank calls
  initMdpOptions() for the modules, and assigns the values from the tree to the
  Options object.  After this, the modules will be exactly in the same state as
  in grompp.

* When other tools (gmx dump or gmx check in particular) read the tpr file,
  they read the structured tree.  In principle, they could operate directly on
  this tree (and `gmx dump` in particular does, with the `-orgir` option).
  However, in the future with proper tpr backward compatibility support, they
  need to call to the modules to ensure that the tree has the structure that
  this version expects, instead of what the original version that wrote the
  file had.  Currently, these tools only call initMdpOptions() and do the basic
  default+native value conversion.

* Any code that is not interested in the parameters for these modules can just
  read the t_inputrec from the tpr file and ignore the tree.

* For compatibility with old tpr files that did not yet have the structured
  tree, the I/O code converts old values for the modules to parameters in the
  structured tree (in tpxio.cpp).

Currently, there is no mechanism for changing the mdp input parameters (adding
new or removing old ones) that would maintain tpr and mdp backward
compatibility.  The vision for this is to use the same transformation engine as
for initMdpTransform() to support specifying version-to-version conversions for
any changed options, and applying the necessary conversions in sequence.  The
main challenge is keeping track of the versions to know which conversions to
apply.

Callbacks to modules during setup and simulation
------------------------------------------------

During setup and simulation, modules receive required information like topologies
and local atom sets by subscribing to callback functions.

To include a notification for your module

* Add the function signature for the callback function to the
  `MDModulesNotifiers` in `mdmodulesnotifiers.h`,

  ```C++
    BuildMDModulesNotifier<...,
                      YourCallbackSignature,
                      ...,
  ```

  (keep alphabetical order for ease of git merge)

* Hand the notifier_ member of the MDModules Implementation class to your
  builder createYourModule(&notifier_)

* Add the function you want to subscribe with in the builder,
  `notifier->subscribe(yourFunction)`

  * To subscribe class member functions of your module, you can use lambda expressions

  ```C++
    notifier->notifier_.subscribe([modulePointer = yourModule.get()]
      (YourCallbackSignature argument){modulePointer(argument);});
  ```

* During setup in , e.g., within `Mdrunner` use

  ```C++
    YourCallbackSignature argument();
    mdModules_.notifier().notifier_.notify(argument);
  ```

Storing non-mdp option module parameters
----------------------------------------

Some mdrun modules want to store data that is non-mdp input, e.g., the result of
computation during setup. Atom indices of index groups are one example:
they are evaluated from strings during grompp time and stored as list of
integers in the run input file. During the mdrun setup the information to
evaluate the index groups is no longer available.

To store parameters, subscribe to the `KeyValueTreeBuilder*` notification that
provides a handle to a KeyValueTreeBuilder that allows adding own information to
that tree.

To restore parameters, subscribe to the `const KeyValueTreeObject &`
notification that returns the tree that is build by the KeyValueTreeBuilder*.
