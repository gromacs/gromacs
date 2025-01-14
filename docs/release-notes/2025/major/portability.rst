Portability
^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Added support to compile |Gromacs| using AMD HIP as GPU backend
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

It is now possible to use AMD HIP directly as the GPU backend to run
simulations on AMD devices. For now only the NBNxM kernels are can
be offloaded to the device using this backend.

:issue:`4947`

Added support for the oneMath interface library for GPU FFTs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This enables cross-vendor support for GPU FFTs to the |Gromacs|
SYCL backend. Either cuFFT or rocFFT can now be used with
Intel DPC++ and Codeplay's plugins for NVIDIA and AMD GPUs.

:issue:`4744`

Update of required build tool versions
""""""""""""""""""""""""""""""""""""""

Updated minimal required versions:

 - CMake 3.28,
 - Compilers: GCC 11, Clang 14 (using Clang 18 for code formatting),
 - GPU toolkits: CUDA 12.1, oneAPI 2024.0, AdaptiveCpp 23.10, ROCm 5.2.

:issue:`5014`
