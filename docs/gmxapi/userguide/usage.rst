========================
Using the Python package
========================

After installing GROMACS, sourcing the "GMXRC" (see GROMACS docs), and installing
the gmxapi Python package (see :doc:`install`), import the package in a Python
script or interactive interpreter. This documentation assumes a convenient alias
of ``gmx`` to refer to the ``gmxapi`` Python package.

::

    import gmxapi as gmx

For full documentation of the Python-level interface and API, use the ``pydoc``
command line tool or the :py:func:`help` interactive Python function, or refer to
the :doc:`pythonreference`.

Python does not wrap any command-line tool, so once installation is complete,
there shouldn't be any additional configuration necessary, and any errors that
occur should be caught at the Python level. Any Python *exception* raised by gmxapi
should be descended from (and catchable as) :class:`gmxapi.exceptions.Error`.

As an exception to the preceding paragraph, we have a tool specifically for
wrapping arbitrary (unintegrated) command line tools: See :class:`gmxapi.commandline_operation`.
