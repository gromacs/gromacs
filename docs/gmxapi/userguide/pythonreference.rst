==============================
gmxapi Python module reference
==============================

.. Concise reference documentation extracted directly from code.
.. For new and non-backwards-compatible features, API versions must be given.

The Gromacs Python package includes a high-level scripting interface implemented
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
