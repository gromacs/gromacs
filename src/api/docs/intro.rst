=====
About
=====

Goals and scope
===============

A (set of) well-designed API(s), compartmentalized from the core Gromacs code,
provide external interfaces to support a modern scientific software ecosystem.

* Flexible workflow
* Flexible compute resources
* Flexible tool selection
* Streamlined code development
* Sustainable code sharing

Design and implementation can and should be encapsulated into iterative feature enhancements to allow incremental forward progress, but components are unified in an architecture designed to meet a variety of client needs.

* Data manipulation
  * file i/o parsing, editing and conversion
  * editing of data structures (topology, structure, input, trajectory)
* Data analysis
  * third party interaction with Gromacs data without reinventing wheels
  * use of existing Gromacs analysis tools
* Driving simulation workflows
  * scripting high-level tasks
  * "bottling" research methods into tools
* Optimizing nontrivial workflows, pipelines, and tool chains
  * data graph execution
  * optimization based on data locality
  * optimization based on compute locality
  * connect simulation and analysis in a pipeline
* Sustainable interfaces to external software
  * plug-in interface
  * Python as a common language for API implementations
  * High-level C++ API on which to build Gromacs tools and external interfaces.

Layers
======

Features are implemented at one or more distinct interface layers, both to
serve different clients and to support sustainable software engineering.

* high-level Python interface

  * convenience API
  * command-line alternative
  * powerful scripting
  * interoperability with external tools
  * uses API to separate high-level tasks from core code

* middleware

  * simulation and analysis chains / pipelines
  * data graph execution
  * execution environment abstraction

* lower level external access from Python and C++ for

  * more control
  * rapid prototyping
  * extensibility
  * compartmentalization from core library

* simulation plugins, call-backs, and integration of external code into data graph

  * developing new methods
  * flexible selection of the right tool for the job
  * benefit from encapsulated third-party tools

For more details, see :doc:`layers`.

Components
==========

Features can be organized into several categories that (probably) constitute
distinct sections of API.

* File-backed data interfaces
* Microstate data
* Runner and modules
* Execution context abstraction
* Simulation internals

For more detail, refer to :doc:`components`
