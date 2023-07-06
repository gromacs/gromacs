Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

HeFFTe multi-GPU FFT plan options are now configurable
""""""""""""""""""""""""""""""""""""""""""""""""""""""

New environment variables ``GMX_HEFFTE_RESHAPE_ALGORITHM``,
``GMX_HEFFTE_USE_GPU_AWARE``, ``GMX_HEFFTE_USE_PENCILS``, and
``GMX_HEFFTE_USE_REORDER`` permit the HeFFTe plan options to be
configured at run time. The performance obtained can vary with the
quality of implementation of e.g. the GPU-aware MPI library, as well
as the layout and number of the GPUs participating in the 3D-FFT.
Users can now find and use the best settings for their case. See
the HeFFTe documentation for more details.
