#include "pybind11/pybind11.h"

/// Test for non-copy of buffer data
/* Make sure that when we expose a buffer to Python and pass it back in,
the raw data has the same address. */

/// Test for correct reference counting (don't destroy buffer providers)
/* Make sure that when a buffer is returned to Python by an object who's
Python references are abandoned, the C++ reference count is maintained
for the life of the buffer. To be extra thorough, monitor from a shared_ptr
copy in another C++ object and not the original shared_ptr returned through the
bindings. Do we need keep_alive? Does it do anything with buffers? */
