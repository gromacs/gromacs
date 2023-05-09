Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

gmxapi.commandline_operation environment variable filtering
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A new utility (:py:func:`gmxapi.runtime.filtered_mpi_environ()`) is available
to remove MPI-related environment variables from :py:data:`os.environ`, such as
to prepare the subprocess environment of `gmxapi.commandline_operation`.

This is a follow-up to :issue:`4423`, for which the original fix appeared to be insufficient.

:issue:`4736`
