Removed functionality
^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!


Support for version 1 of the hardware locality library ``hwloc`` is removed
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Building with ``-DGMX_HWLOC=ON`` now requires ``hwloc`` 2.0 or later,
which has been supported in |Gromacs| for several years.

