Portability
^^^^^^^^^^^

Python environment
""""""""""""""""""

Where Python is required,
`CPython <https://www.python.org>`__ versions 3.6 to 3.8 are supported.

CMake now detects Python using
`FindPython3 <https://cmake.org/cmake/help/v3.13/module/FindPython3.html>`__.
If you previously used ``PYTHON_EXECUTABLE`` to hint the location of the Python
interpreter, you should instead specify the Python "root" or "prefix" path
(the directory containing ``./bin/python3``) with CMake variable
``Python3_ROOT_DIR`` or ``CMAKE_PREFIX_PATH``. As other infrastructure evolves,
``PYTHON_EXECUTABLE`` may cease to have the desired effect without warning.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

CMake
"""""

Updated required CMake version to 3.13.

C++ standard
""""""""""""

|Gromacs| has updated the required C++ standards compliance from C++14 to C++17,
and requires 2017 standard library features. See the install guide for details.

Cygwin
""""""

|Gromacs| now builds on Cygwin with both gcc and clang compilers.

Windows
"""""""

|Gromacs| now builds correctly on Windows with MSVC even when the path
to the source or build directory has a space in it.

Builds with MSVC 2019 correctly detect the proper static linking setup
during CMake configuration.

RDTSCP usage and reporting
""""""""""""""""""""""""""

|Gromacs| now defaults always on x86 to use the RDTSCP machine
instruction for lower latency timing. Very old machines might need to
configure with ``GMX_USE_RDTSCP=off``. Non-x86 platforms are
unaffected, except that they will no longer report that RDTSCP is
disabled (because that is self-evident).

armv8+sve support (ARM_SVE)
"""""""""""""""""""""""""""
Support for ARM Scalable Vector Extensions (SVE) has been added.
|Gromacs| supports SVE vector length fixed at CMake configure time
(typically via the -msve-vector-bits=<len> compiler option),
which is at the time of the release supported in GNU GCC 10 and later,
and will supported soon by LLVM 12 and compilers based on this.
The default is to detect the default vector length at CMake configure time,
and that can be changed with the ``GMX_SIMD_ARM_SVE_LENGTH=<bits>`` option.
Supported values are 128, 256, 512 and 1024. Note that the nonbonded
kernels have not been optimized for ARM_SVE as of yet.
ARM_SVE support is contributed by the Research Organization for Science Information and Technology (RIST)
