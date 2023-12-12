.. _density-guided-simulation:

Applying forces from three-dimensional densities
------------------------------------------------

In density-guided simulations, additional forces are applied to atoms that depend
on the gradient of similarity between a simulated density and a reference density.

By applying these forces protein structures can be made to "fit" densities
from, e.g., cryo electron-microscopy. The implemented approach extends the ones
described in \ :ref:`192 <refOrzechowski2008>`, and \ :ref:`193 <refIgaev2019>`.

Overview
^^^^^^^^

The forces that are applied depend on:

 * The forward model that describes how atom positions are translated into a
   simulated density, :math:`\rho^{\mathrm{sim}}\!(\mathbf{r})`.
 * The similarity measure that describes how close the simulated density is to
   the reference density, :math:`\rho^{\mathrm{ref}}`, :math:`S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\mathbf{r})]`.
 * The scaling of these forces by a force constant, :math:`k`.

The resulting effective potential energy is

.. math:: U = U_{\mathrm{forcefield}}(\mathbf{r}) - k S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\mathbf{r})]\,\mathrm{.}
          :label: eqndensone

The corresponding density based forces that are added during the simulation are

.. math:: \mathbf{F}_{\mathrm{density}} = k \nabla_{\mathbf{r}} S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\mathbf{r})]\,\mathrm{.}
          :label: eqndenstwo

This derivative decomposes into a similarity measure derivative and a simulated
density model derivative, summed over all density voxels :math:`\mathbf{v}`

.. math:: \mathbf{F}_{\mathrm{density}} = k \sum_{\mathbf{v}}\partial_{\rho_{\mathbf{v}}^{\mathrm{sim}}} S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] \cdot \nabla_{\mathbf{r}} \rho_{\mathbf{v}}^{\mathrm{sim}}\!(\mathbf{r})\,\mathrm{.}
          :label: eqndensthree

Thus density-guided simulation force calculations are based on computing a
simulated density and its derivative with respect to the atom positions, as
well as a density-density derivative between the simulated and the reference
density.

Usage
^^^^^

Density-guided simulations are controlled by setting ``.mdp`` options and
providing a reference density map as a file additional to the ``.tpr``.

All options that are related to density-guided simulations are prefixed with
``density-guided-simulation``.

Setting ``density-guided-simulation-active = yes`` will trigger density-guided
simulations with default parameters that will cause atoms to move into the
reference density.

The simulated density and its force contribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Atoms are spread onto the regular three-dimensional lattice of the reference
density. For spreading the atoms onto the grid, the discrete Gauss transform is
used. The simulated density from atoms at positions :math:`\mathbf{r_i}` at a
voxel with coordinates :math:`\mathbf{v}` is

.. math:: \rho_{\mathbf{v}} = \sum_i A_i \frac{1}{\sqrt{2\pi}^3\sigma^3} \exp[-\frac{(\mathbf{r_i}-\mathbf{v})^2}{2 \sigma^2}]\,\mathrm{.}
          :label: eqndensfour

Where :math:`A_i` is an amplitude that is determined per atom type and may be
the atom mass, partial charge, or unity for all atoms.

The width of the Gaussian spreading function is determined by :math:`\sigma`.
It is not recommended to use a spreading width that is smaller than the
grid spacing of the reference density.

The factor for the density force is then

.. math:: \nabla_{r} \rho_{\mathbf{v}}^{\mathrm{sim}}\!(\mathbf{r}) = \sum_{i} - A_i \frac{(\mathbf{r_i}-\mathbf{v})}{\sigma} \frac{1}{\sqrt{2\pi}^3\sigma^3} \exp[-\frac{(\mathbf{r_i}-\mathbf{v})^2}{2 \sigma^2}]\,\mathrm{.}
          :label: eqndensfive

The density similarity measure and its force contribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are multiple valid similarity measures between the reference density and
the simulated density, each motivated by the experimental source of the
reference density data. For the density-guided simulations in |Gromacs|, the following
measures are provided:

The inner product of the simulated density,

.. math:: S_{\mathrm{inner-product}}[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] =
                \frac{1}{N_\mathrm{voxel}}\sum_{v=1}^{N_\mathrm{voxel}} \rho^{\mathrm{ref}}_v \rho^{\mathrm{sim}}_v\,\mathrm{.}
        :label: eqndenssix

The negative relative entropy between two densities,

.. math:: S_{\mathrm{relative-entropy}}[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] =
           \sum_{v=1, \rho^{\mathrm{ref}}>0, \rho^{\mathrm{sim}}>0}^{N_\mathrm{voxel}} \rho^\mathrm{ref} [\log(\rho^\mathrm{sim}_v)-\log(\rho^\mathrm{ref}_v)]\,\mathrm{.}
        :label: eqndensseven

The cross correlation between two densities,

.. math:: S_{\mathrm{cross-correlation}}[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}] =
           \frac{\sum_{v}\left((\rho_v^{\mathrm{ref}} - \bar{\rho}^{\mathrm{ref}})(\rho_v^{\mathrm{sim}} - \bar{\rho}^{\mathrm{sim}})\right)}
           {\sqrt{\sum_v(\rho_v^{\mathrm{ref}} - \bar{\rho}^{\mathrm{ref}})^2 \sum_v(\rho_v^{\mathrm{sim}} - \bar{\rho}^{\mathrm{sim}})^2}}\mathrm{.}
        :label: eqndenscrosscorr
     

Declaring regions to fit
^^^^^^^^^^^^^^^^^^^^^^^^

A subset of atoms may be chosen when pre-processing the simulation to which the
density-guided simulation forces are applied. Only these atoms generate the
simulated density that is compared to the reference density.

Performance
^^^^^^^^^^^

The following factors affect the performance of density-guided simulations

 * Number of atoms in the density-guided simulation group, :math:`N_{\mathrm{atoms}}`.
 * Spreading range in multiples of Gaussian width, :math:`N_{\mathrm{\sigma}}`.
 * The ratio of spreading width to the input density grid spacing, :math:`r_{\mathrm{\sigma}}`.
 * The number of voxels of the input density, :math:`N_{\mathrm{voxel}}`.
 * Frequency of force calculations, :math:`N_{\mathrm{force}}`.
 * The communication cost when using multiple ranks, that is reflected in a constant :math:`c_{\mathrm{comm}}`.

The overall cost of the density-guided simulation is approximately proportional to

.. math:: \frac{1}{N_{\mathrm{force}}} \left[N_{\mathrm{atoms}}\left(N_{\mathrm{\sigma}}r_{\mathrm{\sigma}}\right)^3 + c_{\mathrm{comm}}N_{\mathrm{voxel}}\right]\,\mathrm{.}
          :label: eqndenseight

Applying force every N-th step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cost of applying forces every integration step is reduced when applying the
density-guided simulation forces only every :math:`N` steps. The applied force
is scaled by :math:`N` to approximate the same effective Hamiltonian as when
applying the forces every step, while maintaining time-reversibility and energy
conservation. Note that for this setting, the energy output frequency must be a
multiple of :math:`N`.

The maximal time-step should not exceed the fastest oscillation period of any
atom within the map potential divided by :math:`\pi`. This oscillation period
depends on the choice of reference density, the similarity measure and the force
constant and is thus hard to estimate directly. It has been observed to be
in the order of picoseconds for typical cryo electron-microscopy data, resulting
in a `density-guided-simulation-nst` setting in the order of 100.

Combining density-guided simulations with pressure coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that the contribution of forces from density-guided simulations to the
system virial are not accounted for. The size of the effect on the
pressure-coupling algorithm grows with the total summed density-guided simulation
force, as well as the angular momentum introduced by forces from density-guided
simulations. To minimize this effect, align your structure to the density before
running a pressure-coupled simulation.

Additionally, applying force every N-th steps does not work with the current
implementation of infrequent evaluation of pressure coupling and the constraint
virial.

Periodic boundary condition treatment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Of all periodic images only the one closest to the center of the density map
is considered.

The reference density map format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reference input for the densities are given in mrc format according to the
"EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014".
Closely related formats like ``ccp4`` and ``map`` might work.

Be aware that different visualization software handles map formats differently.
During simulations, reference densities are interpreted as visualised by ``VMD``.

Output
^^^^^^

The energy output file will contain an additional "Density-fitting" term.
This is the energy that is added to the system from the density-guided simulations.
The lower the energy, the higher the similarity between simulated and reference
density.

Adaptive force constant scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To enable a steady increase in similarity between reference and simulated
density while using as little force as possible, adaptive force scaling
decreases the force constant when similarity increases and vice versa. To avoid
large fluctuations in the force constant, change in similarity is measured
with an exponential moving average that smoothens the time series of similarity
measures with a time constant :math:`tau` that is given in ps. If the exponential
moving average similarity increases, the force constant is scaled down by
dividing by :math:`1+\delta t_{\mathrm{density}}/tau`, where
:math:`\delta t_{\mathrm{density}}` is the time between density guided simulation steps.
Conversely, if similarity between reference and simulated density is decreasing,
the force constant is increased by multiplying by :math:`1+2\delta t_{\mathrm{density}}/tau`.
Note that adaptive force scaling does not conserve energy and will ultimately lead to very high
forces when similarity cannot be increased further.

Mapping input structure to density data with affine transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To align input structure and density data, a transformation matrix
:math:`\mathbf{A}` and shift vector :math:`\mathbf{v_{\mathrm{shift}}}` may be
defined that transform the input structure atom coordinates before evaluating
density-guided-simulation energies and forces, so that 

.. math:: U = U_{\mathrm{forcefield}}(\mathbf{r}) - k S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\mathbf{A} \mathbf{r}+\mathbf{v}_{\mathrm{shift}})]\,\mathrm{.}
          :label: eqndensnine

.. math:: \mathbf{F}_{\mathrm{density}} = k \nabla_{\mathbf{r}} S[\rho^{\mathrm{ref}},\rho^{\mathrm{sim}}\!(\mathbf{A} \mathbf{r}+\mathbf{v}_{\mathrm{shift}})]\,\mathrm{.}
          :label: eqndensten

Affine transformations may be used, amongst other things, to perform

 * rotations, e.g., around the z-axis by an angle :math:`\theta` by using :math:`A=\begin{pmatrix} \cos \theta & -\sin \theta & 0 \\ \sin \theta & \cos \theta & 0 \\ 0 & 0 & 1 \end{pmatrix}`.
 * projection, e.g., onto the z-axis by using :math:`A=\begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}`. This allows density-guided simulations to be steered by a density-profile along this axis.
 * scaling the structure against the density by a factor :math:`s` by using :math:`A=\begin{pmatrix} s & 0 & 0 \\ 0 & s & 0 \\ 0 & 0 & s \end{pmatrix}`. This proves useful when, e.g., voxel-sizes in cryo-EM densities have to be adjusted.
 * and arbitrary combinations of these by matrix multiplications (note that matrix multiplications are not commutative).

Future developments
^^^^^^^^^^^^^^^^^^^

Further similarity measures might be added in the future, along with different
methods to determine atom amplitudes. More automation in choosing a force constant
as well as alignment of the input density map to the structure might be provided.
