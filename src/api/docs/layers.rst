=====================
Implementation layers
=====================

Most functionality is available with different levels of granularity and abstraction, or requires lower level facilities for implementation. API level separation is described here. For documentation organized by functionality, see :doc:`components`.

The API is defined with a C++ implementation and a Python implementation (via additional bindings code).

For performance and abstraction, it is important that API layers are able to communicate
with light-weight objects or references. Data and implementation-specific operators are conveyed
with proxy objects. Transfer of actual data across API layers must be explicit
to clarify data ownership and opportunities for optimization. Tasks and operations
requested by API calls are not necessarily executed until specifically requested
or required by later API
calls.

The Python API implementation (bindings to a thin wrapper of the C++ implementation) is used to build a user-friendly and idiomatic Python interface (module). It's functionality is not necessarily mirrored in the CLI tools or a C++ convenience API.

Additional development in the C++ library will be necessary to cleanly support the public API for external software.

Use cases and scenarios for distinct layers are described in the following
documents, with links to relevant components and feature development.

.. toctree::

    layers/python
    layers/bindings
    layers/middleware
    layers/cppapi
    layers/library
