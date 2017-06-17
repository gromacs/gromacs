Allowed language features
=========================

|Gromacs| uses C99 for C files and C++11 for C++ files. 
C++ has a lot of features, but to keep the source code maintainable and easy to read, 
we will avoid using some of them in Gromacs code. The basic principle is to keep things 
as simple as possible.
For compatiblity certain work-arounds are required because not all compilers support 
these standards fully.

* MSVC supports only a subset of C99 and work-arounds are required in those cases.
* Before 7.0 (partial support in 6.5) CUDA didn't support C++11. Therefore any
  header file which is needed (or likely will be nedded) by CUDA should not use C++11.
* We should be able to use virtually all C++ features outside of the header files
  required by CUDA code (and OpenCL kernels), since we have gradually moved to
  compilers that have full support for C++11.

.. TODO: Copy important points from http://www.gromacs.org/index.php?title=Developer_Zone/Programming_Guide/Allowed_C%2B%2B_Features

C++ Standard Library
--------------------

|Gromacs| code must support the lowest common denominator of C++11 standard library
features available on supported platforms.
Some modern features are useful enough to warrant back-porting.
Consistent and forward-compatible headers are provided in ``src/gromacs/compat/``
as described in the `Library documentation <../doxygen/html-lib/group__group__compatibility.xhtml>`_

