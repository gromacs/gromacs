=================
gmxapi data model
=================

stub

Basic data types, containers
============================

Handles and Futures
===================

Proxies and managed resources
=============================

Operations, factories, and data flow: declaration, definition, and initialization
=================================================================================

..  todo:: Reference https://gitlab.com/gromacs/gromacs/-/issues/2993

Expressing inputs and outputs
-----------------------------

Notes on data compatibility
===========================

Avoid dependencies
------------------

The same C++ symbol can have different bindings in each extension module, so
don't rely on C++ typing through bindings. Need schema for PyCapsules.

Adding gmxapi compatible Python bindings should not require dependency on gmxapi
Python package. Compatibility through interfaces instead of inheritance.

