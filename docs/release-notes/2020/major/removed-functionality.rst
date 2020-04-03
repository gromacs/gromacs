Removed functionality
^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Group cut-off scheme
""""""""""""""""""""

The group cut-off scheme has been removed. Several kinds of simulation
that depend on it no longer work.

   * Simulations under vacuum conditions are not supported.
   * User supplied tables for short-range nonbonded interactions are not supported.
   * Switched short-range nonbonded interactions with PME are not supported. 
   * Membrane embedding is deactivated.
   * QMMM is not supported.

:issue:`1852`

Generalized reaction-field
""""""""""""""""""""""""""

This only worked correctly with the group scheme. Note that generalized
reaction-field simulations can still be performed using standard
reaction field and computing the dielectric constant manually.
       
gmx anadock
"""""""""""
The gmx anadock tool was removed since it does not belong in gromacs
(it analyzes AutoDock outputs).

gmx dyndom
""""""""""
The gmx dyndom tool was removed since it does not belong in gromacs
(it analyzes DynDom outputs).

gmx morph
"""""""""
The gmx morph tool was removed since it yields non-physical structures
that can easily be done by a script.

gmx mdrun -gcom
"""""""""""""""

This feature sometimes overrode the effects of various .mdp settings
in a way that was difficult to understand and report. A user who wants
to do communication between PP ranks less often should choose their
``nst*`` mdp options accordingly.
