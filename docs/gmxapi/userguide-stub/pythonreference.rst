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

.. seealso::

    This copy of the documentation was built without
    :doc:`installing the gmxapi package <../userguide/install>`,
    and therefore lacks the full module reference.
    Refer to :ref:`gmxapi_package_documentation` for instructions on building
    complete documentation, or `view online <http://manual.gromacs.org/current/gmxapi/>`__.
