==============================
gmxapi Python module reference
==============================

.. contents:: :local:
    :depth: 2

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
   :members:
   :member-order: groupwise

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

.. autoexception:: gmxapi._gmxapi.Exception

    Root exception for the C++ extension module. Derives from `gmxapi.exceptions.Error`.

.. autoexception:: gmxapi._gmxapi.NotImplementedError

.. autoexception:: gmxapi._gmxapi.ProtocolError

.. autoexception:: gmxapi._gmxapi.UnknownException

.. autoexception:: gmxapi._gmxapi.UsageError

Functions
---------

Tools for launching simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gmxapi._gmxapi.from_tpr

Tools to manipulate TPR input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gmxapi._gmxapi.copy_tprfile

.. autofunction:: gmxapi._gmxapi.read_tprfile

.. autofunction:: gmxapi._gmxapi.write_tprfile

.. autofunction:: gmxapi._gmxapi.rewrite_tprfile

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
