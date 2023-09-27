# h5xx — a template-based C++ wrapper for the HDF5 library

**Build status**: [![build status](https://git.imp.fu-berlin.de/h5md/h5xx/badges/master/build.svg)](https://git.imp.fu-berlin.de/h5md/h5xx/commits/master)

## Introduction

h5xx is a template-based C++ wrapper for the HDF5 library.  The goal of h5xx is
to provide an easy-to-use yet flexible and powerful interface to HDF5 for C++
codes.  In some sense, h5xx aims at providing similar functionality to C++ as
[h5py](http://www.h5py.org/ "HDF5 for Python") does to Python.  For example, a
NumPy-like slicing notation wrapping HDF5 hyperslabs is implemented.
Currently, h5xx supports std::vector, boost::array, and boost::multi_array
containers.

Using h5xx, a hello world example to write a 32x32 Boost multi_array to a HDF5
dataset labeled "my_array" would look as follows:
```
boost_array_2d_t array(boost::extents[32][32]);
// fill array with data
h5xx::file file("data.h5", h5xx::file::out);
h5xx::create_dataset(file, "my_array", array);
h5xx::write_dataset(file, "my_array", array);
```

The unit tests and example codes provide further information on and guidance to
the usage and the capabilities of h5xx.

h5xx can be used, copied, modified, and distributed freely under the terms of
the 3-clause BSD license.


## Requirements

* A compiler that supports at least C++11 is required.  h5xx was developed
and tested using g++ on x86_64 Linux.  h5xx was tested in addition with Intel
icpc.
* HDF5: h5xx requires an installation of the HDF5 library.  HDF5 may be built
with MPI to support parallel IO.  The (deprecated) HDF5 C++ bindings are *not* required.
* Boost: h5xx requires an installation of the Boost C++ library.
h5xx supports the Boost array and the Boost multidimensional array datatype.
Moreover, h5xx uses the Boost `enable_if` set of templates to control the creation
of SFINAE (substitution-failure-is-not-an-error) conditions.
* CMAKE: The cmake build system is required to compile the examples and unit
tests that are packaged with h5xx.

See the file INSTALL.md for some hints on installing and using h5xx.


## Supported platforms

In 2016/02, h5xx was successfully built and tested on the following
configurations.

* Ubuntu Linux 14.04 (x86_64)
  * gcc: 4.9.3
  * OpenMPI: 1.8.5
  * Boost: 1.{49, 55, 58, 60}
  * HDF5: 1.8.{12, 14, 15, 16}
* SUSE LINUX ENTERPRISE SERVER (SLES) 11 (x86_64)
  * gcc: 4.8.4, 4.9.3
  * Intel Compiler (icpc): 14.0, 15.0
  * Intel MPI: 4.1.3, 5.0.3
  * Boost: 1.57, 1.58
  * HDF5: 1.8.12, 1.8.16

h5xx is considered mature.  It is under active development.
