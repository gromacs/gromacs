New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


The AWH exponential histogram growth can now be controlled
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The accelerated weight histogram growth factor during the initial phase
was hard-coded to 3. Now this value can be controlled by the user.
It is set to 2 by default for increased stability.

If the TPR was generated with an earlier |Gromacs| version,
the old default value of 3 will be used.


Added support for instrumentation based on wallcycle regions using NVTX/ROCTX/ITT
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Basic support has been added for GPU tracing libraries so wallcycle main and sub-regions
will show up in tracing timelines which can help with performance analysis.
The tracing instrumentation support can be enabled with one of the following CMake variables:
``GMX_USE_NVTX``, ``GMX_USE_ROCTX``, ``GMX_USE_ITT``.

:issue:`4446`
