Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Instant-submission mode enabled by default when building with AdaptiveCpp
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In |Gromacs| 2024, one had to manually enable the instant-submission mode
of AdaptiveCpp when building GROMACS
(``-DSYCL_CXX_FLAGS_EXTRA=-DHIPSYCL_ALLOW_INSTANT_SUBMISSION=1``).
Now it is enabled by default, improving performance by up to 20%
when running on GPU and slightly reducing the CPU usage when using
SYCL/AdaptiveCpp backend.

More OpenMP parallelization
"""""""""""""""""""""""""""

OpenMP multithreading is now used for position restraints, GPU LINCS initialization,
domain decomposition state vector sorting and a few other places. This is particularly
relevant for runs where a fast GPU is bottlenecked by DD and other CPU tasks.

:issue:`5100`

GPU buffer operations enabled by default
""""""""""""""""""""""""""""""""""""""""

This is generally faster on modern GPUs, but, for very small systems, this can worsen performance by a few percent.
In this case, the old behavior can be restored by setting ``GMX_GPU_DISABLE_BUFFER_OPS=1`` environment variable.
