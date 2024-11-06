Deprecated functionality
------------------------

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Changes anticipated to |Gromacs| 2025 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| 2025
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MTTK pressure coupling is deprecated
""""""""""""""""""""""""""""""""""""

:issue:`5072`

The TNG trajectory format is deprecated
"""""""""""""""""""""""""""""""""""""""

The TNG file format will be removed in a future release. It will be replaced by the more
widely used HDF5-based H5MD format. There will be at least one version of |Gromacs|
supporting both H5MD and TNG to allow conversions between the two formats.

:issue:`5225`
