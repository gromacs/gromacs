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
