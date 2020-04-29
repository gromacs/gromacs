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
``Python3_ROOT`` or ``CMAKE_PREFIX_PATH``. As other infrastructure evolves,
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

GROMACS has updated the required C++ standards compliance from C++14 to C++17,
and requires 2017 standard library features. See the install guide for details.
