Improvements to |Gromacs| tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

``gmx msd`` has been migrated to the trajectoryanalysis framework
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The tool now uses the |Gromacs| selection syntax. Rather than piping selections via stdin,
selections are now made using the "-sel" option. There is a new option called ``-maxtau``,
which limits maximum time delta between frames to compare for calculating MSDs. This will allow
users who otherwise would run into out-of-memory errors and slow execution with large systems
to restrict sampling to useful tau values.

This migration comes with about a 20% speedup in execution time.

Some rarely used features have yet to be migrated, including:

- The -tensor option is not yet implemented.
- System COM removal with -rmcomm has not yet been implemented.
- B-factor writing using the -pdb option is not yet supported.

A slight behavior change is the removal of the -mw option. ``gmx msd`` with ``-mol`` will
take the MSD of the center-of-mass of of molecules, while no mass-weighting is done
when ``-mol`` is not selected. In previous |Gromacs| versions, ``-mw`` was on by default,
and ``-nomw`` was silently ignored when ``-mol`` was chosen. This change will only cause
different results when performing MSD calculations on a non-homogenous group of particles without
``-mol`` set.

:issue:`2368`

``gmx chi`` no longer needs ``residuetypes.dat`` entries for custom residues
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The need to add the names of custom residues to ``residuetypes.dat`` has been
removed, because it served no purpose. This makes ``gmx chi`` easier to use.

``gmx wham`` has had minor improvements to its text output
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Reporting about file handling and input-file column contents are easier to
follow.
