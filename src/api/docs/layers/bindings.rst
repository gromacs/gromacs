====================
Low-level Python API
====================

Gromacs public API features are available from both a C++ API and a Python API.
The Python API is implemented as a C++ extension, providing a low-level interface
as a submodule in the Gromacs Python package. It is sensible to implement the
high-level Python interface in pure Python because the flexibility of the language
makes it easy to provide a stable interface. For the low-level Python module,
implementation entirely in C++ makes it easy to maintain consistent semantics
for the different language implementations with a minimally cumbersome development
process. To provide an idiomatic Python interface, though, the C++ API is not
(necessarily) directly wrapped with a Python bindings exporter, but some
additional C++ shims are necessary.

The Python API allows convenience APIs, user interfaces, and higher-level tools
to be implemented entirely in Python, using core Gromacs functionality. It
facilitates interoperability with external software that also provides a Python
API.

In addition to supporting efficient interaction, the API allows clients to
receive maximum benefit from the parallelism and performance optimizations in
the core Gromacs library. It provides a stable interface and abstraction level
to decouple Gromacs code from external projects.

Performance and compatibility can best be guaranteed for API calls operating
exclusively on API objects. Gromacs input, output, data structures, and objects
are managed as API objects. Objects defined on one side of an API call are
mapped to references or proxy objects on the other side, transfering data only
when necessary or explicitly requested. To facilitate data locality optimization,
computation is deferred across a sequence of API calls, when possible.

The API defers not only computation and data transfer, but also control flow
related to execution context. This allows task scheduling and parallelism to
be managed separately, either at a low level by the library or at a higher level
through the middleware layer.
.. seealso:: middleware
Though agnostic to the execution context, the API provides client-side tools
with which to interact with the execution context where necessary, such as to
support data locality optimizations in client code.

Bindings requirements
=====================

The bindings introduce no additional external dependencies for build, installation,
run-time linking, or Python use. The module supports, but doesn't require,
data sharing with tools like ``numpy`` through the Python buffer protocol.

For ease of maintenance and to help assure consistent behavior in Python and C++
API, the exported Python interface, the bindings tools, and any shims necessary to
adapt the C++ API are all implemented in C++ at the standards level required by
Gromacs.

Any bindings tools must be compatible with the Gromacs license.

The bindings manage passing and sharing object references back and forth between
the Python interpreter and the C++ API with propoer reference counting and
lifetime management.

Suggestion
==========

`Pybind11 <http://pybind11.readthedocs.io/en/stable/intro.html>`_

It is designed specifically to interface between Python and C++11 (or
greater). Only C++ syntax is needed to implement bindings, making it much
simpler and lighter-weight than more general tools for interfacing between more
than a pair of languages.

It is small and efficient, taking full advantage of
the Python C API and modern C++ capabilities to allow the Python interpreter to
serve as the glue to connect C++ objects, providing an alternative to dynamic
loading from C++.

It uses syntax and idioms familiar to anyone who previously used Boost Python,
but is much lighter weight: it does not link against any other libraries or
have any additional dependencies at build time or run time.
