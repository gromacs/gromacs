Miscellaneous
^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

grompp no longer modifies nstcomm
"""""""""""""""""""""""""""""""""

grompp will no longer set nstcomm, the interval for center of mass motion
removal, equal to nstcalcenergy when nstcomm < nstcalcenergy.
A note is still printed in that case.

Bonded atom types names can now start with a digit
""""""""""""""""""""""""""""""""""""""""""""""""""

Bonded atom types names in topologies were not allowed to start with a number.
Now all names are supported that contain at least one non-digit character.

:issue:`4120`

grompp now warns when exclusion forces might be missing
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

When using PME, exclusions between non-perturbed,atom pairs should be within
the cut-off distance, otherwise mdrun might not compute grid correction forces
and energies. grompp now computes these distance for the starting structure
and warns when they are beyond 90% of the cut-off distance and generates
an error when they are beyond the cut-off distance.

:issue:`4051`

Core spin-up code is removed
""""""""""""""""""""""""""""""""""""""""""""""""""

Formerly, on non-x86 and non-PowerPC platforms, mdrun ran some
multi-threaded code to try to wake up any cores that the OS might have
powered down. This caused problems on some Arm platforms, and does not
seem to suit a significant number of platforms for use of GROMACS. So
now it is removed.

If required, please manually spin-up the cores with, e.g., ``stress --cpu $(nproc --all)``.

:issue:`4074`

Add documentation for linear angle potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Added documentation and reference for the linear angle potential. Also
added please_cite entry, but there is no call to reference it yet.

:issue:`4286`

gmxapi.mdrun guarantees trajectory output
"""""""""""""""""""""""""""""""""""""""""

gmxapi simulations now always run with full-precision trajectory output (``-o``)
in order to guarantee the availability of a usable output trajectory through the
``mdrun.output.trajectory`` result.

:issue:`4285`

gmxapi.mdrun accepts arbitrary runtime arguments
""""""""""""""""""""""""""""""""""""""""""""""""

Arbitrary mdrun arguments can be passed through gmxapi with the new *runtime_args* key word
argument, accepting a dictionary of flags and values.

:issue:`4284`
