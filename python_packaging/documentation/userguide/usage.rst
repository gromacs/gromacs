========================
Using the Python package
========================


For full documentation of the Python-level interface and API, use the ``pydoc``
command line tool or the ``help()`` interactive Python function, or refer to
the :ref:`python-procedural` documentation.

Python does not wrap any command-line tool, so once installation is complete,
there shouldn't be any additional configuration necessary, and any errors that
occur should be caught at the Python level. Exceptions should all be descendants
of :class:`gmx.exceptions.Error`.

Running simulations
===================

Plugins
-------

The work graph
--------------

Parallelism
-----------

Synchronous and asynchronous work
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Accelerated computation
^^^^^^^^^^^^^^^^^^^^^^^

Simulation input
================

Preparing simulation input
--------------------------

With work graph operations
^^^^^^^^^^^^^^^^^^^^^^^^^^

With file utilities
^^^^^^^^^^^^^^^^^^^

Manipulating simulation parameters
----------------------------------

With key word arguments
^^^^^^^^^^^^^^^^^^^^^^^

With work graph operations
^^^^^^^^^^^^^^^^^^^^^^^^^^

Accessing trajectory data
=========================

From files
----------

From in-memory objects
----------------------

From gmxapi operation results
-----------------------------
