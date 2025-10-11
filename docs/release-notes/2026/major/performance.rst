Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

Optimized PME force transfers for multi-GPU runs with GPU-aware MPI
===================================================================

When using a GPU build with an (external) GPU-aware MPI library, the transfer of PME force data between
GPUs has been optimized to use non-blocking MPI calls, improving performance by up to around 5%.


.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

More permissive AWH histogram equilibration check
"""""""""""""""""""""""""""""""""""""""""""""""""

The AWH histogram equilibration check, which is recommended to be used with
multiple walkers, used a too strict tolerance of 20% maximum deviation.
This tolerance can now be set by the user and the default value is 50%.
This value provide significantly faster convergence when using multiple walkers.
Since it does not deteriorate the convergence with a single walkers,
the histogram equilibration check is now turned on by default.
