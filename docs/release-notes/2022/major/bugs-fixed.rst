Bugs fixed
^^^^^^^^^^

Fixed slight inaccuracies when using virtual sites with pressure coupling
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Virtual sites were reconstructed after the system was propagated, but before
scaling due to pressure coupling. For virtual site types which are not a linear
combination of other atoms, this is not completely correct. Since the scaling
due to pressure coupling is very small in healthy simulations, the resulting
inaccuracies are expected to have been extremely minor, and in most cases
undetectable.

:issue:`3866`

Correct dVremain/dl when nstdhdl > nstcalcenergy
""""""""""""""""""""""""""""""""""""""""""""""""

When nstcalcenergy was not a multiple of nstdhdl, incorrect dVremain/dl
terms were written in the energy file. Note that all dH/dl output in
both dhdl.xvg and the energy file, which is used by e.g. gmx bar, was correct.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Removed velocity output for acceleration groups
"""""""""""""""""""""""""""""""""""""""""""""""

The reported velocity in the energy file for acceleration groups was always
zero. Now their velocity is no longer reported in the energy file.

:issue:`1354`

Use correct c0 parameter in Me2PO4 in OPLSAA
""""""""""""""""""""""""""""""""""""""""""""

OPLSAA torsions must sum to 0, but the parameters for Me2PO4 did not do so. Changed the c0
parameter to the correct value.

:issue:`4075`

Allow function type Fourier Dihedral with free energy perturbations
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Fourier Dihedral (dihedral interaction type 3) could not be used in
free energy perturbation simulations. Under the hood the dihedral parameters
were anyhow converted to Ryckaert-Bellemans parameters, so now the checks
for perturbations are the same for the two functions.

:issue:`2606`

Do not scale coordinates of frozen atoms during Parrinello-Rahman pressure coupling
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When Parrinello-Rahman pressure coupling was used, the box scaling was applied to all the atoms,
causing frozen atoms to shift. The effect is more drastic towards the sides of the box and when the
pressure is changed significantly during the simulations. Now, the frozen atoms will be ignored by
the coupling and atoms with frozen dimensions shall keep such values.

:issue:`3075`

Avoid non-uniform rotation with Test Particle Insertion in anisotropic systems
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With anisotropic systems the random angles would not get a uniform distribution.

:issue:`3558`

Allow free energy calculations with a linear angle potential
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Free energy calculations with a linear angle potential were not
explicitly allowed by grompp.

:issue:`3456`


Fixed progress display in trjconv and trjcat
""""""""""""""""""""""""""""""""""""""""""""

The progress information (frame number and time) shown during trajectory 
operations in trjconv and trjcat is now correctly displayed.

:issue:`4320`

Fixed GROMOS dihedral generation for disulfide bridges
""""""""""""""""""""""""""""""""""""""""""""""""""""""

The pdb2gmx functionality now generates correct dihedrals for disulfide
bridges with the GROMOS force field series.

:issue:`4188`

Fixed energy term naming for periodic improper dihedrals
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Those used the same name internally as the non-periodic version for printing
to energy files and reading from them. This could cause tools being confused
when trying to compare terms from files where the terms where written in
a different order.

gmx density now always uses relative coordinates
""""""""""""""""""""""""""""""""""""""""""""""""

There is no realistic use case for using absolute coordinates in binning
when the box dimension is changing, so gmx density now always uses
relative coordinates internally. This also avoids issues with output
scaling to the last instead of average box size when users forget
this option, ensures the output is always correct, and gets rid of
occassional segfaults.

:issue:`3830`
