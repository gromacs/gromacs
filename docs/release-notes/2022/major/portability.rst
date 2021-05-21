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
