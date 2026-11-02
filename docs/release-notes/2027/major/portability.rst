Portability
^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

clFFT is no longer accepted on macOS
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Apple OpenCL compiler fails to produce a functional clFFT build, so
configuring an OpenCL build with ``-DGMX_GPU_FFT_LIBRARY=clFFT`` on macOS
is now rejected by CMake, rather than leaving PME on the GPU unavailable
at run time. VkFFT, which is the default on macOS, should be used instead.

Libtorch versions that require C++20 work with the default C++17 build
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Libtorch headers from PyTorch 2.14 onwards require C++20, which used to
stop the build unless all of |Gromacs| was configured with
``-DCMAKE_CXX_STANDARD=20``. Now only the sources that include them (the
NNPot interface and Colvars) are compiled with C++20; the rest of |Gromacs|
keeps the configured standard.

:issue:`5754`
