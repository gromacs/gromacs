========================
Python packaging roadmap
========================

Describe functional goals, feature tracking, and associated tests.

Functionality to build
======================

General tool implementation sequence
------------------------------------

1. Tools in ``operation`` submodule wrap Python code with compliant interface
2. ``commandline_operation()`` provides UI for (nearly) arbitrary CLI tools
3. ``tool`` submodule provides trivially wrapped ``gmx`` tools with well-defined inputs and outputs

Simulation tool implementation sequence
---------------------------------------

1. ``mdrun``
2. ``modify_input``
3. ``make_input``

Data flow
---------

* ``result()`` methods force data locality for use by non-API calls
* output properties allow proxied access to (fault tolerant, portable) result futures
* ``logical_and()``, etc. allow manipulation of data / signals "in flight"

Data flow topology tools
------------------------

* implicit scatter / map
* implicit broadcast
* ``gather()`` makes results local
* ``reduce()`` allows operation across an ensemble (implicit allReduce where appropriate)

Control flow
------------

* ``subgraph`` allows several operations to be bundled with a scope that allows
  user-declared data persistence (e.g. when used in a loop)
* ``while_loop`` allows repeated execution of a subgraph, subject to a declared
  condition of the subgraph state/output, managing accessibilty of data handles
  internal and external to the graph

Execution environment
---------------------

1. Top level module function ``run()`` launches work to produce requested data
2. Infrastructure manages / negotiates dispatching or allocation of available resources

Packaging
---------

1. Python package can be installed from GROMACS source directory after GROMACS installation.
2. Python package can be installed by a user from a hand-built GROMACS installation.
3. Python package can be installed from a Linux binary GROMACS distribution using
   appropriate optimized binary libraries.

Installation procedures are testable with Dockerfiles in a Docker-based CI environment.
Ref: https://github.com/kassonlab/gmxapi/blob/master/ci_scripts/test_installers.sh

Typing
------

In addition to basic scalar types,
some structured data types are immediately necessary.

* gmxapi.NDArray supports an N-dimensional array of a single scalar gmxapi data type
  (with sufficient metadata for common array view APIs)
* gmxapi.Map supports an associative container with key strings and Variant values

NDArray is immediately necessary to disambiguate non-scalar parameters/inputs from
implicit data flow topological operations.

Data structures currently coupled to file formats should be decomposable into
Maps with well-specified schema, but in the absence of more complete API data
abstractions we just use filename Strings as the API-level content of
input/output data handles. (Suggest: no need for file objects at
the higher-level API if read_Xfiletype(/write_Yfiletype) operations
produce(/consume) named and typed data objects of specified API types.)

The API specification should be clear about policies for narrowing and widening
of scalar types.

Sequence
========

Features development sequence based on functional priorities and dependencies.

* fr1: wrap importable Python code.
* fr2: output proxy establishes execution dependency (superseded by fr3)
* fr3: output proxy can be used as input
* fr4: dimensionality and typing of named data causes generation of correct work topologies
* fr5: explicit many-to-one or many-to-many data flow
* fr7: Python bindings for launching simulations
* fr8: gmx.mdrun understands ensemble work
* fr9: MD plugins
* fr10: fused operations for use in looping constructs
* fr11: Python access to TPR file contents
* fr12: Simulation checkpoint handling
* fr13: ``run`` module function simplifies user experience
* fr14: Easy access to GROMACS run time parameters
* fr15: Simulation input modification
* fr16: Create simulation input from simulation output
* fr17: Prepare simulation input from multiple sources
* fr18: GROMACS CLI tools receive improved Python-level support over generic commandline_operations
* fr19: GROMACS CLI tools receive improved C++-level support over generic commandline_operations
* fr20: Python bindings use C++ API for expressing user interface
* fr21 User insulated from filesystem paths
* fr22 MPI-based ensemble management from Python
* fr23 Ensemble simulations can themselves use MPI

Expectations on Mark for Q1-Q2 2019 GROMACS master changes
==========================================================

* Broker and implement build system amenable to multiple use
  cases. Need to be able to build and deploy python module from single
  source repo that is usable (i.e. can run the acceptance tests).

  - Some kind of nested structure likely appropriate, perhaps
    structured as nested CMake projects that in principle could stand
    alone. That's probably workable because nested projects can see
    the parent project's cache variables (TODO check this)
  - probably a top-level project coordinating a libgromacs build and a
    python module build, with the former typically feeding the latter
  - the libgromacs build may be able to leverage independent efforts
    towards a multi-configuration build (so SIMD/MPI/GPU agnostic)
  - top-level project offers much the same UI as now, passing much of
    it through to the libgromacs project
  - top-level project offers the option to find a Python (or be told
    which to use), to find a libgromacs (or be told, or be told to
    build), to build any necessary wrapper binaries (ie. classical gmx
    and mdrun), and to deploy all linked artefacts to
    CMAKE_INSTALL_PREFIX or the appropriate Python site-packages
  - the top-level project will be used by e.g. setup.py wrapper
    from scikit-build/distutils
  - requires reform of compiler flags handling
  - probably requires some re-organization of external dependencies
    of libgromacs
  - follow online "Modern CMake" best practices as far as practicable
  - library should be available for static linking with position
    independent code to allow a single shared object to be built for
    the Python module.

* Dissolve boundary between libgmxapi and libgromacs

  - no effort on form and stability of the C++ headers and library in
    2019, beyond what facilitates implementing the Python interface
    in GROMACS 2020
  - existing libgromacs declarations of "public API" and installed
    headers removed

* libgromacs to be able to be use an MPI communicator passed in,
  rather than hard-coding MPI_COMM_WORLD anywhere. It is likely that
  existing wrapper binaries can use the same mechanism to pass
  MPI_COMM_WORLD to libgromacs.

* UI helpers should express.
  - preferred name for datum as a string: ``nsteps``, ``tau-t``, etc.
  - setter (function object, pointer to a builder method, )
  - typing and type discovery (could be deducible from setter, but something to allow user input checking, or determination
    of the suitability of a data source to provide the given input)
  - help text: can be recycled to provide auto-extracted documentation, command-line help, and annotation in Python docstrings.
  - for CLI: short name for flag. E.g. 'p' for "topology_file"
  - for compatibility: deprecated / alternate names. E.g. "nstlist" for "neighbor_list_rebuild_interval", or "orire" for
    "enable_orientation_restraints"
  - default values

Possible GROMACS source changes whose impact is currently unknown
=================================================================
* gmx::Any (which is a flavour of C++17 std::any) type could be
  helpful at API boundary. Also perhaps a flavour of C++17
  std::optional or std::variant.

Additional goals
================

Some project goals are integrations or optimizations that are explicitly hidden from the user
and not testable in a high level script, but should be reflected as milestones in a roadmap.

GROMACS source changes deferred to later in 2019
================================================
* Build system works also from tarball
* Build system can produce maximally static artefacts (for performance
  on HPC infrastructure)
* express grompp and mdrun options handling with gmx::Options to
  prepare for future dictionary-like handling in Python without
  serializing a .tpr file
