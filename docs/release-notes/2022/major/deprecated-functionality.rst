Deprecated functionality
------------------------

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Changes anticipated to |Gromacs| 2022 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| 2022
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GMX_OPENCL_NB_CLUSTER_SIZE CMake variable deprecated in favor of GMX_GPU_NB_CLUSTER_SIZE
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Both OpenCL and SYCL support different cluster sizes, so GMX_GPU_NB_CLUSTER_SIZE should
be used going forward.

The built-in viewer ``gmx view`` will be removed
""""""""""""""""""""""""""""""""""""""""""""""""

There is little use and no tests of this functionality, so it is not worth attempting to
maintain moving forward.

:issue:`4296`

The analysis tool ``gmx chi`` will be removed
"""""""""""""""""""""""""""""""""""""""""""""

This tool has not been functional for a few years.
Please comment at the linked issue if you have any interest in it.

:issue:`4108`

Guessing masses and atomic radii from atom names is deprecated
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When atom masses or van-der-Waals radii are needed, we suggest building
a proper |Gromacs| topology instead of using PDB files directly, even
if the tool supports it.

:issue:`3368`
:issue:`4288`

