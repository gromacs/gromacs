==============================
gmxapi Python module reference
==============================

.. Concise reference documentation extracted directly from code.
.. For new and non-backwards-compatible features, API versions must be given.

The |Gromacs| Python package includes a high-level scripting interface implemented
in pure Python and a lower-level API implemented as a C++ extension module.
The pure Python implementation provides the basic ``gmxapi`` module and
classes with a very stable syntax that can be maintained with maximal compatibility
while mapping to lower level interfaces that may take a while to sort out. The
separation also serves as a reminder that different execution contexts may be
implemented quite diffently, though Python scripts using only the high-level
interface should execute on all.

Package documentation is extracted from the ``gmxapi`` Python module and is also available
directly, using either ``pydoc`` from the command line or :py:func:`help` from within Python, such
as during an interactive session.

Refer to the Python source code itself for additional clarification.

.. seealso:: :ref:`gmxapi_package_documentation`

.. Configuration for doctest: automated syntax checking in documentation code snippets
.. testsetup::

    import gmxapi as gmx
    from gmxapi.data import tpr_filename

Interface concepts
==================

.. py:currentmodule:: gmxapi.abc

*gmxapi* commands return *references* to operations.
Generally, the operations are collected into a graph of data flow dependencies,
and only executed when the results are requested.

.. autoclass:: OperationReference
    :members:

*gmxapi* uses a `Future` to reference an operation output or data that may not yet be available.

.. autoclass:: Future
    :members:

An `OperationReference` may provide several named Futures on its *output* attribute.

A `Future` may be provided directly as inputs to other *gmxapi* commands.
*gmxapi* will execute the required operation to get the data when it is needed.

To get an actual result in your Python script, you can call
:py:func:`~Future.result()` on any *gmxapi* data reference.
If the operation has not yet been executed, the operation (and any operation dependencies)
will be executed immediately.

You can also force an operation to run by calling its :py:func:`~OperationReference.run()` method.
But this is not generally necessary unless your only goal is to produce output files on disk
that are not consumed in the same script.

In some cases, a `Future` can be subscripted to get a new Future representing
a slice of the original.
For instance, `commandline_operation` objects have a *file* output that produces
a mapping of command line flags to output files (per the *output_files* parameter).
This *file* output can be subscripted with a single command line option to get
a `Future` for just one output file type.
See :ref:`gmxapi simulation preparation` for an illustrative example.

Ensemble data flow
------------------

*gmxapi* automatically generates arrays of operations and parallel data flow,
when parallel inputs are provided to *gmxapi* command parameters.

When a `Future` represents the output of an ensemble operation,
:py:func:`~Future.result()` returns a list with elements corresponding
to the ensemble members.

It is not currently possible to get a `Future` for a specific ensemble member.

See :ref:`gmxapi ensemble` for more information.

gmxapi basic package
=====================

::

    import gmxapi as gmx

.. automodule:: gmxapi

.. autodecorator:: function_wrapper

.. autofunction:: commandline_operation

.. autofunction:: subgraph

.. autofunction:: while_loop

Simulation module
=================

.. automodule:: gmxapi.simulation

.. py:currentmodule:: gmxapi

Preparing simulations
---------------------

.. autofunction:: read_tpr

.. For the OutputDataProxy classes, the creation signatures are distracting.
   As of Sphinx 6.0, however, the autodoc_class_signature configuration option
   cannot be overridden for individual 'autoclass' directives, and in general
   the default "mixed" value seems appropriate for this documentation.

.. autoclass:: gmxapi.simulation.read_tpr.OutputDataProxy

.. autofunction:: modify_input

.. autoclass:: gmxapi.simulation.modify_input.OutputDataProxy

Running simulations
-------------------

.. autofunction:: mdrun

.. autoclass:: gmxapi.simulation.mdrun.OutputDataProxy

Utilities
=========

.. automodule:: gmxapi.utility

.. autofunction:: config

.. autofunction:: join_path

.. py:currentmodule:: gmxapi

.. autofunction:: concatenate_lists

.. autofunction:: join_arrays

.. autofunction:: logical_not

.. autofunction:: make_constant

Status messages and Logging
===========================

.. automodule:: gmxapi._logging
   :members:

Exceptions module
=================
..  automodule:: gmxapi.exceptions
    :members:

gmx.version module
==================
..  automodule:: gmxapi.version
    :members:

Core API
========

.. automodule:: gmxapi._gmxapi

Exceptions
----------

Module Exceptions
~~~~~~~~~~~~~~~~~

.. autoexception:: gmxapi._gmxapi.Exception

    Root exception for the C++ extension module. Derives from `gmxapi.exceptions.Error`.

.. autoexception:: FeatureNotAvailable


Wrapped C++ exceptions emitted through the supporting |Gromacs| library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoexception:: gmxapi._gmxapi.MissingImplementationError

.. autoexception:: gmxapi._gmxapi.ProtocolError

.. autoexception:: gmxapi._gmxapi.UsageError

Other
~~~~~

No other C++ exceptions are expected, but will be wrapped in a
:py:class:`Exception` to help tracing and reporting bugs.

.. autoexception:: gmxapi._gmxapi.UnknownException

Functions
---------

This documentation is provided for completeness and as an aid to developers.
Users of the :py:mod:`gmxapi` package, generally, should not need to use the
following tools directly.

Tools for launching simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gmxapi._gmxapi.from_tpr

.. autofunction:: gmxapi._gmxapi.create_context

Tools to manipulate TPR input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gmxapi._gmxapi.copy_tprfile

.. autofunction:: gmxapi._gmxapi.read_tprfile

.. autofunction:: gmxapi._gmxapi.write_tprfile

.. autofunction:: gmxapi._gmxapi.rewrite_tprfile

Utilities
~~~~~~~~~

.. autofunction:: gmxapi._gmxapi.has_feature

    Available features may depend on the package version, the details of the
    supporting |Gromacs| installation, the software environment detected
    when the package was built, or possibly on detected runtime details.
    These feature checks are largely for internal use. The :py:mod:`gmxapi`
    commands may adjust their behavior slightly depending on feature checks,
    and (at worst) should produce meaningful error messages or exceptions.

    Named features:

    * *create_context*: `create_context` can be used to initialize a `Context`
      with assigned resources.
    * *mpi_bindings*: C++ extension module was built with :py:mod:`mpi4py` compatibility.


Classes
-------

.. autoclass:: gmxapi._gmxapi.Context
    :members:

.. autoclass:: gmxapi._gmxapi.MDArgs
    :members:

.. autoclass:: gmxapi._gmxapi.MDSession
    :members:

.. autoclass:: gmxapi._gmxapi.MDSystem
    :members:

.. autoclass:: gmxapi._gmxapi.SimulationParameters
    :members:

.. autoclass:: gmxapi._gmxapi.TprFile
    :members:
