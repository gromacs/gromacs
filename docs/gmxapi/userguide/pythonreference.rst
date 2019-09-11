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

The following documentation is extracted from the ``gmxapi`` Python module and is also available
directly, using either ``pydoc`` from the command line or :py:func:`help` from within Python, such
as during an interactive session.

Refer to the Python source code itself for additional clarification.

.. Configuration for doctest: automated syntax checking in documentation code snippets
.. testsetup::

    import gmxapi as gmx
    from gmxapi.data import tpr_filename

.. _python-procedural:

gmxapi basic pacakage
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

Python API
==========

.. contents:: :local:

Python Classes
--------------

Python context managers
-----------------------

gmx.exceptions module
---------------------
..  automodule:: gmxapi.exceptions
    :members:

gmx.version module
------------------
..  automodule:: gmxapi.version
    :members:

Core API
========

.. automodule:: gmxapi._gmxapi

Functions
---------

Classes
-------
