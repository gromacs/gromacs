Shear simulations
-----------------

A common type of non-equilibrium simulations in fluid dynamics and rheology are
shearing simulations. These are non-equilibrium simulations where work is
performed on the simulation system to achieve a shear flow. This can be used
to compute viscosities and friction and to study the effect of shear stress on conformations.
In |Gromacs| there are four different ways to achieve shear flow.

Groups of atoms can be given a constant acceleration, which is effectively
a mass-weighted force. This will cause such groups to move with respect to
the rest of the system. Care needs to be taken to control the velocity of
the center of mass of the system. Normal center of mass motion removal
can not be used, as that would affect the flow in the system.

As |Gromacs| supports general triclinic unit-cell shapes, the unit cell can
be deformed to set up a shear flow. This can be achieved either by deforming
the unit cell directly using the ``deform`` option in the :ref:`mdp` file,
or this can be driven by applying an off-diagonal stress through pressure
coupling. In the former case, one can measure the viscosity through
the stress, in the latter case through measuring the shear rate.

For measuring the viscosity of simple liquids one can use a cosine-shaped
acceleration profile, which can be specified using the ``cos-acceleration``
option in the :ref:`mdp` file. As the unit-cell does not deform, this
avoids some complications of the other methods. The viscosity is computed
on the fly and reported in the energy file.

And finally, there is the case where one wants to study the effect of walls
on the flow. In particular, structured walls are of interest, consisting
of atoms that can be of any kind. In this case one wants to have walls
on two sides of the system, typically in the xy-plane close to z=0 and
the box height. The flow is then driven by moving the walls at constant
speed by using a constant force. A constant force can be achieved by
use of acceleration groups, but that will not allow position restraining
atoms in the walls along the direction of the shear, which is needed
for some types of walls. For the case of walls where (part of) the atoms
are position restrained, a constant speed can be set by using
the free-energy lambda-coupling code. To achieve this, you need to supply
a second, B-state, position restraint file with the ``-r`` option
of :ref:`gmx grompp`. If you shift the coordinates in this file by 1 nm
in the direction of shear, you can set the speed of the walls with the
``delta-lambda`` option in the :ref:`mdp` file. Note that this makes
lambda increase proportionally with simulation time. There is no limit
on magnitude of lambda and periodic shifts of walls are handled
correctly. When the position restraint coordinates are shifted by 1 nm,
the force on the walls is given directly by :math:`dV/d\lambda`.

A Poiseuille flow is a popular setup in experiments. Unfortunately this is
difficult to achieve in simulations. The best would be to, as in experiment,
apply a pressure difference over (part of) the simulation box. But that
is not easy to set up. One can accelerate all liquid atoms, but this does
not guarantee that atoms that interact directly with the wall experience
the same forces as they would in an experiment. A slightly better setup
would be accelerating only the atoms in the middle of the flow,
but spatially defined acceleration groups are currently not supported
in |Gromacs|.
