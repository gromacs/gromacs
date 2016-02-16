Allowed language features
=========================

GROMACS uses C99 for C files and C++11 for C++ files. 
C++ has a lot of features, but to keep the source code maintainable and easy to read, 
we will avoid using some of them in Gromacs code. The basic principle is to keep things 
as simple as possible.
For compatiblity certain work-arounds are required because not all compilers support 
these standards fully.

* MSVC supports only a subset of C99 and work-arounds are required in those cases.
* Before 7.0 (partial support in 6.5) CUDA didn't support C++11. Therefore any
  header file which is needed (or likly will be nedded) by CUDA should not use C++11.
* C++11 features which are not widely implemented (including in MSVC 2015 and GCC 4.6)
  should not be used.

.. TODO: Copy important points from http://www.gromacs.org/index.php?title=Developer_Zone/Programming_Guide/Allowed_C%2B%2B_Features
