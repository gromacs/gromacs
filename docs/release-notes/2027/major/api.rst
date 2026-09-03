Changes to the API
^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Direct storage of ``TimeControl`` affects functions in ``pargs.h`` and ``trxio.h``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

During command line argument parsing, time control flags (``-b``, ``-e`` and ``-dt``)
are set to an input ``gmx::TimeControl`` struct (in ``gromacs/fileio/timecontrol.h``)
instead of to global variables. The following functions are affected:

``parse_common_args`` takes an extra ``gmx::TimeControl*`` input, which can be ``nullptr``.
If time control is desired a non-``nullptr`` object must be passed.

``read_first_frame`` and ``read_first_x`` both take an extra ``gmx::TimeControl*`` input,
which can be ``nullptr`` if time control is not desired.

``check_times`` takes an extra ``gmx::TimeControl&`` input containing the time control
parameters to check against.
