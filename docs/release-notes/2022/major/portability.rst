Portability
^^^^^^^^^^^

Intel classic compiler (icc/icpc) no longer supported
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

We now support the Intel clang-based compiler from oneAPI (icx/icpx)
instead. Please use it, or gcc.

:issue:`3893`

Provisional: Initialize GMX_INSTALL_NBLIB_API and GMXAPI build options from BUILD_SHARED_LIBS
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

CMake options ``GMXAPI`` and ``GMX_INSTALL_NBLIB_API`` produce shared object libraries,
so their default values are now initialized from ``BUILD_SHARED_LIBS``.
Pending movement on :issue:`3605` and related issues, the coupling between these
options is subject to change, but users generally should not need to manually set
``GMXAPI`` and ``GMX_INSTALL_NBLIB_API``.

:issue:`4053`

Updates to pybind11 dependency
""""""""""""""""""""""""""""""

pybind11 is no longer bundled with |Gromacs|.

The gmxapi 0.3 Python package build system relies on PEP 517/518 build requirements to get pybind11 header
dependencies through the Python packaging system. Package managers like ``pip`` will download dependencies
automatically. Package managers that do not automatically fulfill dependencies should still report the missing
dependency to the user.

The ``sample_restraint`` sample project
(bundled in ``python_packaging/sample_restraint``)
still has a primitive CMake-only build procedure.
If you fork a project from this source, you may choose to modernize the build system (similarly to that of
``gmxapi``) or to bundle the pybind11 sources.
Within the GROMACS repository, the ``sample_restraint`` option default is now ``GMXAPI_EXTENSION_DOWNLOAD_PYBIND=ON``.

:issue:`4092`

Bundle muparser
"""""""""""""""

|Gromacs| now bundles MuParser version 2.3. It is also possible
to link to an external provided library.
