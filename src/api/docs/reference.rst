=========
Reference
=========

Concise reference documentation extracted directly from code.
For more curated instruction, see the :doc:`userguide` and :doc:`developerguide`.

.. For new and non-backwards-compatible features, API versions must be given.

Python module reference
=======================

documentation generated with ``sphinx.ext.doctest``

The Gromacs Python interface is implemented as a high-level scripting interface implemented in pure Python and a lower-level API implemented as a C++ extension.
The pure Python implementation provides the basic ``gmx`` module and
classes with a very stable syntax that can be maintained with maximal compatibility
while mapping to lower level interfaces that may take a while to sort out. The
separation also serves as a reminder that different execution contexts may be
implemented quite diffently, though Python scripts using only the high-level
interface should execute on all. Bindings to the ``libgromacs`` C++ API are
provided in the submodule ``gmx.core``.

C++ API reference
=================

Documentation extracted for public C++ API. All symbols and names appearing in
the public API headers must be documented to the extent of their visibility to
external code. Design details are described in the :doc:`developerguide`
and more implementation details may be available separately in the library
documentation. If the API components are built in developer mode, links are
provided for classes that provide additional interfaces to library code. All
documentation must indicate relevant API versions.
