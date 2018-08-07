This is a fork of the [main Gromacs project](http://www.gromacs.org/) in which interface, API, and extensibility issues are being investigated.
The forked project lives on GitHub at [https://github.com/kassonlab/gromacs-gmxapi](https://github.com/kassonlab/gromacs-gmxapi/)

[![Build Status](https://travis-ci.org/kassonlab/gromacs-gmxapi.svg?branch=master)](https://travis-ci.org/kassonlab/gromacs-gmxapi)

In addition to a regular GROMACS installation, this fork provides `libgmxapi` for
high-level C++ access to GROMACS MD simulation.
It exists primarily to support the [`gmxapi`](https://github.com/kassonlab/gmxapi) companion project that provides a Python module and bindings.

This README.md file supplants the main README file to avoid merge conflicts while providing convenient documentation to the repository browser.

# Installation

Install as you would a regular copy of GROMACS. The following example downloads the source into a directory named `gromacs`,
creates a parallel (out-of-source) `build` directory, configures, builds, and installs. Use e.g. `make -j10 install` to build in parallel with 10 processes.

    $ git clone https://github.com/kassonlab/gromacs-gmxapi.git
    $ mkdir build
    $ cd build
    $ cmake ../gromacs-gmxapi -DCMAKE_INSTALL_PREFIX=/path/to/where/i/want/gromacs -DGMX_THREAD_MPI=ON -DGMX_GPU=OFF -DGMX_BUILD_OWN_FFTW=ON                                                                                                                                                                                             
    $ make install

You may then either source the GMXRC file (as usual for GROMACS use) or export the environment variable
`gmxapi_DIR=/path/to/where/i/want/gromacs` to help `gmxapi` clients such as the Python 
package or your own CMake project to find
what it needs to build against the gmxapi library.

# Documentation

To build additional documentation, use the additional build targets `gmxapi_cppdocs` and `gmxapi_cppdocs_dev`.
You will need to have `doxygen` installed.
If you would like to write code that uses `libgmxapi`, use `make gmxapi_cppdocs`.
For more detail, or if you would like to extend or contribute to the API, `make gmxapi_cppdocs_dev`.

Then refer either to `docs/html/doxygen/api-user/index.html` or
`docs/html/doxygen/api-dev/index.html` in the build directory.

Also, please use the issue tracking system or feel free to suggest other modes of communication.

# Releases, Compatibility, and Versioning

Beginning with the 1.0 release, gmxapi will follow traditional semantic versioning for API compatibility.
ABI compatibility guarantees are TBD. For 0.0.x, expect API incompatibility for each release. For 0.1.x
and onwards, patch releases should maintain API backwards compatibility with the minor version.

The version numbers for gmxapi are encoded in the repository solely in the `src/api/CMakeLists.txt` file.
During CMake configuration, the `src/api/cpp/gmxapi/version.h` file is created so that the built library can
report this version through the `gmxapi::Version` interface. Client code should include the installed 
`gmxapi/version.h` header in order to embed the constants `gmxapi::GMXAPI_MAJOR`, `gmxapi::GMXAPI_MINOR`,
and `gmxapi::GMXAPI_PATCH` so that API compatibility checks can be performed at runtime.

When a new software release is tagged, the next commit on the development branch should increment the patch level to distinguish development builds from the tagged release. As incompatibilities are introduced
in feature branches, minor or major version number should be incremented as appropriate. At this time,
client code has no indication of whether the version presented in a development build of gmxapi is an
officially specified API revision or is subject to change. Developers coding against development branches
should keep this in mind. If this becomes problematic, please offer your suggestions or propose a revision
to the `gmxapi::Version` API.
