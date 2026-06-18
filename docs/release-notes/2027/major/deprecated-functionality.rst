Deprecated functionality
------------------------

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Changes anticipated to |Gromacs| 2027 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| 2027
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Double-batched FFT (BBFFT) backend
""""""""""""""""""""""""""""""""""

The BBFFT backend has been deprecated in |Gromacs| following Intel's abandonment of the project.
Users on Intel GPU hardware are advised to rely on the default oneMKL backend.

:issue:`5658`

