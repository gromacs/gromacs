Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

The deform option was unsuited for flow simulations
"""""""""""""""""""""""""""""""""""""""""""""""""""

The deform option now only deforms the box and does not modify atom positions
anymore. In contrast to previous versions, it instead corrects the velocities
of particles when they are shifted by a periodic box vector. Now, deform is
also useful for shear flows. Applications where the system was stretched until
some interactions broke were probably not affected measurably by
these issues. Note that a velocity profile should be generated when using
deform with the current or later versions. An mdp option has been added
to let ``grompp`` do this.

:issue:`4607`

mdrun now checks for excluded pairs beyond the cut-off with reaction-field and FEP
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With reaction-field electrostatics and free-energy calculations,
excluded atom pairs are not allowed to be beyond the Coulomb cut-off distance.
Now mdrun checks for this and throws an error when this occurs.

:issue:`4667`

Limit pressure deviations due to missing Lennard-Jones interactions
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For systems dominated by Lennard-Jones interactions, i.e. with no or very weak
electrostatics, e.g. most coarse-grained systems, the Verlet buffer was often
set such that missing Lennard-Jones interactions could lead to the pressure
increasing by more than 1 bar over the lifetime of the pair list. Now an mdp
parameter has been added to limit the deviation in the average pressure.
The default tolerance is 0.5 bar.

:issue:`4861`

enemat now prints correct headers when using ``-free`` or ``-eref`` options
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Fixed a long-standing bug when ``gmx enemat`` would output incorrect headers
to XVG.

:issue:`4812`

`gmxapi.commandline_operation` implicitly converts *input_files* to absolute paths
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Relative paths in the *input_files* mapping are now explicitly relative to the working
directory of the caller.

:issue:`4827`
