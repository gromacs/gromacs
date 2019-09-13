.. _density-guided-simulation:

Applying forces from three-dimensional densities
------------------------------------------------

In density-guided simulations, additional forces are applied to atoms that depend
on the gradient of similarity between a simulated density and a reference density.

By applying these forces protein structures can be made to "fit" densities
from, e.g., cryo electron-microscopy. The implemented approach extends the ones 
described in \ :ref:`182 <refOrzechowski2008>`, and \ :ref:`183 <refIgaev2019>`.

Overview
^^^^^^^^

The forces that are applied depend on: 

 * The forward model that describes how atom coordinates are translated into a 
   simulated density, :math:`\rho^{\mathrm{sim}}\!(\vec{r})`. 
   
 * The similarity measure that describes how close the simulated density is to 
   the reference density, :math:`\rho^{\mathrm{ref}}`, :math:`S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\vec{r})]`.

 * The scaling of these forces by a force constant, :math:`k`.

The resulting effective potential energy is 

.. math:: 
    U = U_{\mathrm{forcefield}}(\vec{r}) + k S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\vec{r})]\,\mathrm{.}

The corresponding density based forces that are added during the simulation are

.. math::
    \vec{F}_{\mathrm{density}} = -k \nabla_{\vec{r}} S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\vec{r})]\,\mathrm{.}

This derivative decomposes into a similarity measure derivative and a simulated
density model derivative, summed over all density voxels :math:`\vec{v}`

.. math::
    \vec{F}_{\mathrm{density}} = -k \sum_{\vec{v}}\partial_{\rho_{\vec{v}}^{\mathrm{sim}}} S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] \cdot \nabla_{\vec{r}} \rho_{\vec{v}}^{\mathrm{sim}}\!(\vec{r})\,\mathrm{.}

Thus density-guided simulation force calculations are based on computing a
simulated density and its derivative with respect to the atom coordinates, as
well as a density-density derivative between the simulated and the reference
density.

Usage
^^^^^

Density-guided simulations are controlled by setting `.mdp` options and
providing a reference density map as a file additional to the `.tpr`.

All options that are related to density-guided simulations are prefixed with 
`density-guided-simulation`.

Setting `density-guided-simulation-active = yes` will trigger density-guided
simulations with default parameters that will cause atoms to move into the 
reference density.

The simulated density and its force contribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Atoms are spread onto the regular three-dimensional lattice of the reference
density. For spreading the atoms onto the grid, the discrete Gauss transform is
used. The simulated density from atoms at positions :math:`\vec{r_i}` at a
voxel with coordinates :math:`\vec{v}` is

.. math:: 
    \rho_{\vec{v}} = \sum_i A_i \frac{1}{\sqrt{2\pi}^3\sigma^3} \exp[-\frac{(\vec{r_i}-\vec{v})^2}{2 \sigma^2}]\,\mathrm{.}

Where :math:`A_i` is an amplitude that is determined per atom type and may be
the atom mass, partial charge, or unity for all atoms.

The width of the Gaussian spreading function is determined by :math:`\sigma`.
It is not recommended to use a spreading width that is smaller than the
grid spacing of the reference density.

The factor for the density force is then

.. math::
    \nabla_{r} \rho_{\vec{v}}^{\mathrm{sim}}\!(\vec{r}) = \sum_{i} - A_i \frac{(\vec{r_i}-\vec{v})}{\sigma} \frac{1}{\sqrt{2\pi}^3\sigma^3} \exp[-\frac{(\vec{r_i}-\vec{v})^2}{2 \sigma^2}]\,\mathrm{.}

The density similarity measure and its force contribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are multiple valid similarity measures between the reference density and
the simulated density, each motivated by the experimental source of the
reference density data. For the density-guided simulations in |Gromacs|, the following
measures are provided:

The inner product of the simulated density,

.. math:: S_{\mathrm{inner-product}}[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] =
                \frac{1}{N_\mathrm{voxel}}\sum_{v=1}^{N_\mathrm{voxel}} \rho^{\mathrm{ref}}_v \rho^{\mathrm{sim}}_v\,\mathrm{.}

The negative relative entropy between two densities,

.. math:: S_{\mathrm{relative-entropy}}[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] =
           \sum_{v=1, \rho^{\mathrm{ref}}>0, \rho^{\mathrm{sim}}>0}^{N_\mathrm{voxel}} \rho^\mathrm{ref} [\log(\rho^\mathrm{sim}_v)-\log(\rho^\mathrm{ref}_v)]\,\mathrm{.}

The cross correlation between two densities,

.. math:: S_{\mathrm{cross-correlation}}[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] =
           \frac{\sum_{v}\left((\rho_v^{\mathrm{ref}} - \bar{\rho}^{\mathrm{ref}})(\rho_v^{\mathrm{sim}} - \bar{\rho}^{\mathrm{sim}})\right)}
           {\sqrt{\sum_v(\rho_v^{\mathrm{ref}} - \bar{\rho}^{\mathrm{ref}})^2 \sum_v(\rho_v^{\mathrm{sim}} - \bar{\rho}^{\mathrm{sim}})^2}}\mathrm{.}

     

Declaring regions to fit
^^^^^^^^^^^^^^^^^^^^^^^^

A subset of atoms may be chosen when pre-processing the simulation to which the
density-guided simulation forces are applied. Only these atoms generate the
simulated density that is compared to the reference density.

Applying force every nth step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cost of applying forces every integration step is reduced when applying the
density-guided simulation forces only every :math:`N` steps. The applied force
is scaled by :math:`N` to approximate the same effective Hamiltonian as when 
applying the forces every step, while maintaining time-reversibility and energy
conservation.

The maximal time-step should not exceed the fastest oscillation period of any 
atom within the map potential divided by :math:`\pi`. This oscillation period
depends on the choice of reference density, the similarity measure and the force
constant and is thus hard to estimate directly. It has been observed to be
in the order of picoseconds for typical cryo electron-microscopy data, resulting
in a `density-guided-simulation-nst` setting in the order of 100.

Periodic boundary condition treatment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Of all periodic images only the one closest to the center of the density map
is considered.

The reference density map format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reference input for the densities are given in mrc format according to the 
"EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014".
Closely related formats like `ccp4` and `map` might work.

Be aware that different visualization software handles map formats differently.
During simulations, reference densities are interpreted as visualised by `VMD`.
If the reference map shows unexpected behaviour, swapping endianess with a map
conversion tool like `em2em` might help.

Output
^^^^^^

The energy output file will contain an additional "Density-fitting" term. 
This is the energy that is added to the system from the density-guided simulations.
The lower the energy, the higher the similarity between simulated and reference
density.

Future developments
^^^^^^^^^^^^^^^^^^^

Further similarity measures might be added in the future, along with different
methods to determine atom amplitudes. More automation in choosing a force constant
as well as alignment of the input density map to the structure might be provided.
