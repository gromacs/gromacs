Portability
^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Added support for the oneMKL interface library for GPU FFTs
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This enables cross-vendor support for GPU FFTs to the |Gromacs|
SYCL backend. Either cuFFT or rocFFT can now be used with
Intel DPC++ and Codeplay's plugins for NVIDIA and AMD GPUs.

:issue:`4744`
