New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Density-guided simulations
""""""""""""""""""""""""""

Users can now apply additional forces from three dimensional reference
densities. These forces can be used to "fit" atoms into the densities by
increasing the similarity of a simulated density to the reference density.

Multiple protocols are available for how to calculate simulated densities
as well as how the similarity between a reference and a simulated density is
evaluated.

Virtual site on the line through two atoms at fixed distance
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This is use useful for e.g. halogens in the CHARMM force field.

:issue:`2451`

gmxapi Python support
"""""""""""""""""""""

Data flow driven simulation and analysis from Python is now available in a
default |Gromacs| installation when users install the *gmxapi* Python package.
See :ref:`gmxapi`.

New modular simulator
"""""""""""""""""""""

A new approach for how to combine individual calculation steps during a single simulation
step is introduced, with focus on extensibility and modularization. This simulator is the default
now for simulations using velocity-verlet in NVE, NVT (v-rescale thermostat only), NPT (v-rescale
thermostat and Parrinello-Rahman barostat only), or NPH (Parrinello-Rahman barostat only), with or
without free energy perturbation.
