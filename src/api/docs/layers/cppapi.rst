=======
C++ API
=======

Any CLI tools worth maintaining should be implemented as much as possible in
terms of the public API. It should be straight-forward to implement module runners
as well as modules without dependence on CLI, filesystem access, or terminal I/O.

Reference counted objects should be used so that it is straight-forward for API
clients to extend the life of data objects without necessitating memory copies.
E.g. API clients can only trust gmx::ArrayRef and t_trxframe objects for the
duration of an API call when the client implements a specified interface.

File formats should have public APIs.

ABI and API compatibility policies should be discussed and described. Gromacs
and many other scientific software packages are distributed as source rather
than binary, so the need to recompile plugins after building a new Gromacs
doesn't seem that onerous, nor does requiring a plugin to be built with the
same compiler and settings as extracted from FindGROMACS.cmake, but users who
want to tie together tools from several projects into new code may find
themselves with irreconcilable requirements. It may be sufficient to accept
potential STL compatibilities, say, at the C++ level and rely on Python (or a C
API) to glue together more finicky code with raw pointers and/or copy conversion
of incompatible data classes.
