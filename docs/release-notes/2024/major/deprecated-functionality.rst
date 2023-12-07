Deprecated functionality
------------------------

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Changes anticipated to |Gromacs| 2024 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The analysis tool ``gmx chi`` no longer deprecated
""""""""""""""""""""""""""""""""""""""""""""""""""

Given the community interest, the decision was made to keep ``gmx chi``.

:issue:`4108`

The analysis tool ``gmx gyrate-legacy`` deprecated
""""""""""""""""""""""""""""""""""""""""""""""""""

``gmx gyrate`` has been partly re-implemented in the modern |Gromacs| analysis framework.
The old implementation is still available as ``gmx gyrate-legacy``. Please plan to update to the new version.
Please let the |Gromacs| developers know of desired functionality missing from, or broken in, the new implementation.


:issue:`4927`

The analysis tool ``gmx hbond-legacy`` deprecated
"""""""""""""""""""""""""""""""""""""""""""""""""

``gmx hbond`` has been partly re-implemented in the modern |Gromacs| analysis framework.
The old implementation is still available as ``gmx hbond-legacy``. Please plan to update to the new version.
Please let the |Gromacs| developers know of desired functionality missing from, or broken in, the new implementation.

:issue:`4927`

The analysis tools ``gmx sans`` and ``gmx saxs`` deprecated
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``gmx sans`` and ``gmx saxs`` has been partly re-implemented under new name ``gmx scattering`` in the modern |Gromacs| analysis framework.
The old implementations are still available as ``gmx sans-legacy`` and ``gmx saxs-legacy``. Please plan to update to the new version.
Please let the |Gromacs| developers know of desired functionality missing from, or broken in, the new implementation.

:issue:`4927`

Functionality deprecated in |Gromacs| 2024
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Xeon Phi support will be removed
""""""""""""""""""""""""""""""""""""

Intel Xeon Phi series of accelerators has been discontinued in 2020,
and most supercomputers using it are now retired. The support for
Xeon Phi is deprecated in |Gromacs| 2024 and will be removed in the
next release.

:issue:`4740`

