# Installation of h5xx

h5xx is a C++ template library and does not require compilation. To use h5xx in
your code, include `h5xx.hpp` in your source file and add the h5xx directory to
the compiler's include search path (using `-I`). In addition, make sure that the
HDF5 and Boost headers and the HDF5 libraries are passed to the compiler/linker.

To compile the unit tests and examples, the script `h5xx_build.sh` can be used
as a starting point. CMAKE is required to build the executables.  The CMAKE
configuration files give an example of how to integrate h5xx, HDF5 and Boost
into a code project.

## Troubleshooting

  - *Linking to the Boost C++ library fails.* Check that the paths to the Boost
    C++ headers and the library are consistent, in particular if you are using a
    custom Boost installation. Toggling the CMake flag "Boost_USE_MULTITHREADED"
    and redetecting the Boost C++ installation may resolve the issue.
