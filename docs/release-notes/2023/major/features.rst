New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

mdrun now also reports the conserved energy quantity with AWH bias sharing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

gmxapi.mdrun now publishes the simulation working directory path
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``gmxapi.mdrun().output.directory`` provides the (Future) path to the
working directory/directories for the simulation(s).
This can be useful in conjunction with :py:func:`gmxapi.utility.join_path`
to express data flow based on files known to the user to be produced by
the simulation but not represented by other existing attributes of the
:py:class:`~gmxapi.simulation.mdrun.OutputDataProxy`.

:issue:`4548`
