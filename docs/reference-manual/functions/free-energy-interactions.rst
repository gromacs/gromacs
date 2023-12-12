.. _feia:

Free energy interactions
------------------------

This section describes the :math:`\lambda`-dependence of the potentials
used for free energy calculations (see sec. :ref:`fecalc`). All common
types of potentials and constraints can be interpolated smoothly from
state A (:math:`\lambda=0`) to state B (:math:`\lambda=1`) and vice
versa. All bonded interactions are interpolated by linear interpolation
of the interaction parameters. Non-bonded interactions can be
interpolated linearly or via soft-core interactions.

Starting in |Gromacs| 4.6, :math:`\lambda` is a vector, allowing different
components of the free energy transformation to be carried out at
different rates. Coulomb, Lennard-Jones, bonded, and restraint terms can
all be controlled independently, as described in the
:ref:`mdp` options.

Harmonic potentials
~~~~~~~~~~~~~~~~~~~

The example given here is for the bond potential, which is harmonic in
|Gromacs|. However, these equations apply to the angle potential and the
improper dihedral potential as well.

.. math:: \begin{aligned}
          V_b     &=&{\frac{1}{2}}\left[{(1-{\lambda})}k_b^A + 
                          {\lambda}k_b^B\right] \left[b - {(1-{\lambda})}b_0^A - {\lambda}b_0^B\right]^2  \\
          {\frac{\partial V_b}{\partial {\lambda}}}&=&{\frac{1}{2}}(k_b^B-k_b^A)
                          \left[b - {(1-{\lambda})}b_0^A + {\lambda}b_0^B\right]^2 + 
          		\nonumber\\
                  & & \phantom{{\frac{1}{2}}}(b_0^A-b_0^B) \left[b - {(1-{\lambda})}b_0^A -{\lambda}b_0^B\right]
          		\left[{(1-{\lambda})}k_b^A + {\lambda}k_b^B \right]\end{aligned}
          :label: eqnfepharmpot

GROMOS-96 bonds and angles
~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourth-power bond stretching and cosine-based angle potentials are
interpolated by linear interpolation of the force constant and the
equilibrium position. Formulas are not given here.

Proper dihedrals
~~~~~~~~~~~~~~~~

For the proper dihedrals, the equations are somewhat more complicated:

.. math:: \begin{aligned}
          V_d     &=&\left[{(1-{\lambda})}k_d^A + {\lambda}k_d^B \right]
                  \left( 1+ \cos\left[n_{\phi} \phi - 
          		    {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B
          		    \right]\right)\\
          {\frac{\partial V_d}{\partial {\lambda}}}&=&(k_d^B-k_d^A) 
                   \left( 1+ \cos
          		 \left[
          		    n_{\phi} \phi- {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B
          		 \right]
          	 \right) +
          	 \nonumber\\
                  &&(\phi_s^B - \phi_s^A) \left[{(1-{\lambda})}k_d^A - {\lambda}k_d^B\right] 
                  \sin\left[  n_{\phi}\phi - {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B \right]\end{aligned}
          :label: eqnfeppropdihedral

**Note:** that the multiplicity :math:`n_{\phi}` can not be
parameterized because the function should remain periodic on the
interval :math:`[0,2\pi]`.

Tabulated bonded interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For tabulated bonded interactions only the force constant can
interpolated:

.. math:: \begin{aligned}
                V  &=& ({(1-{\lambda})}k^A + {\lambda}k^B) \, f \\
          {\frac{\partial V}{\partial {\lambda}}} &=& (k^B - k^A) \, f\end{aligned}
          :label: eqnfeptabbonded

Coulomb interaction
~~~~~~~~~~~~~~~~~~~

The Coulomb interaction between two particles of which the charge varies
with :math:`{\lambda}` is:

.. math:: \begin{aligned}
          V_c &=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[{(1-{\lambda})}q_i^A q_j^A + {\lambda}\, q_i^B q_j^B\right] \\
          {\frac{\partial V_c}{\partial {\lambda}}}&=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[- q_i^A q_j^A + q_i^B q_j^B\right]\end{aligned}
          :label: eqnfepcoloumb

where :math:`f = \frac{1}{4\pi \varepsilon_0} = {138.935\,458}` (see
chapter :ref:`defunits`).

Coulomb interaction with reaction field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Coulomb interaction including a reaction field, between two
particles of which the charge varies with :math:`{\lambda}` is:

.. math:: \begin{aligned}
          V_c     &=& f\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]
          \left[{(1-{\lambda})}q_i^A q_j^A + {\lambda}\, q_i^B q_j^B\right] \\
          {\frac{\partial V_c}{\partial {\lambda}}}&=& f\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]
          \left[- q_i^A q_j^A + q_i^B q_j^B\right]
          \end{aligned}
          :label: eqdVcoulombdlambda

**Note** that the constants :math:`k_{rf}` and :math:`c_{rf}` are
defined using the dielectric constant :math:`{\varepsilon_{rf}}` of the
medium (see sec. :ref:`coulrf`).

Lennard-Jones interaction
~~~~~~~~~~~~~~~~~~~~~~~~~

For the Lennard-Jones interaction between two particles of which the
*atom type* varies with :math:`{\lambda}` we can write:

.. math:: \begin{aligned}
          V_{LJ}  &=&     \frac{{(1-{\lambda})}C_{12}^A + {\lambda}\, C_{12}^B}{{r_{ij}}^{12}} -
          \frac{{(1-{\lambda})}C_6^A + {\lambda}\, C_6^B}{{r_{ij}}^6}   \\
          {\frac{\partial V_{LJ}}{\partial {\lambda}}}&=&\frac{C_{12}^B - C_{12}^A}{{r_{ij}}^{12}} -
          \frac{C_6^B - C_6^A}{{r_{ij}}^6}
          \end{aligned}
          :label: eqdVljdlambda

It should be noted that it is also possible to express a pathway from
state A to state B using :math:`\sigma` and :math:`\epsilon` (see
:eq:`eqn. %s <eqnsigeps>`). It may seem to make sense physically to vary the
force field parameters :math:`\sigma` and :math:`\epsilon` rather than
the derived parameters :math:`C_{12}` and :math:`C_{6}`. However, the
difference between the pathways in parameter space is not large, and the
free energy itself does not depend on the pathway, so we use the simple
formulation presented above.

Kinetic Energy
~~~~~~~~~~~~~~

When the mass of a particle changes, there is also a contribution of the
kinetic energy to the free energy (note that we can not write the
momentum :math:`\mathbf{p}` as
:math:`m \mathbf{v}`, since that would result in the
sign of :math:`{\frac{\partial E_k}{\partial {\lambda}}}` being
incorrect \ :ref:`99 <refGunsteren98a>`):

.. math:: \begin{aligned}
          E_k      &=&     {\frac{1}{2}}\frac{\mathbf{p}^2}{{(1-{\lambda})}m^A + {\lambda}m^B}        \\
          {\frac{\partial E_k}{\partial {\lambda}}}&=&    -{\frac{1}{2}}\frac{\mathbf{p}^2(m^B-m^A)}{({(1-{\lambda})}m^A + {\lambda}m^B)^2}\end{aligned}
          :label: eqnfepekin

after taking the derivative, we *can* insert
:math:`\mathbf{p}` = :math:`m \mathbf{v}`, such that:

.. math:: {\frac{\partial E_k}{\partial {\lambda}}}~=~    -{\frac{1}{2}}\mathbf{v}^2(m^B-m^A)
          :label: eqnfepekinderivative

Constraints
~~~~~~~~~~~

The constraints are formally part of the Hamiltonian, and therefore they
give a contribution to the free energy. In |Gromacs| this can be
calculated using the LINCS or the SHAKE algorithm. If we have
:math:`k = 1 \ldots K` constraint equations :math:`g_k` for LINCS, then

.. math:: g_k     =       | \mathbf{r}_{k} | - d_{k}
          :label: eqnfepconstr

where :math:`\mathbf{r}_k` is the displacement vector
between two particles and :math:`d_k` is the constraint distance between
the two particles. We can express the fact that the constraint distance
has a :math:`{\lambda}` dependency by

.. math:: d_k     =       {(1-{\lambda})}d_{k}^A + {\lambda}d_k^B
          :label: eqnfepconstrdistdep

Thus the :math:`{\lambda}`-dependent constraint equation is

.. math:: g_k     =       | \mathbf{r}_{k} | - \left({(1-{\lambda})}d_{k}^A + {\lambda}d_k^B\right).
          :label: eqnfepconstrlambda

The (zero) contribution :math:`G` to the Hamiltonian from the
constraints (using Lagrange multipliers :math:`\lambda_k`, which are
logically distinct from the free-energy :math:`{\lambda}`) is

.. math:: \begin{aligned}
          G           &=&     \sum^K_k \lambda_k g_k    \\
          {\frac{\partial G}{\partial {\lambda}}}    &=&     \frac{\partial G}{\partial d_k} {\frac{\partial d_k}{\partial {\lambda}}} \\
                      &=&     - \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}
          :label: eqnconstrfreeenergy

For SHAKE, the constraint equations are

.. math:: g_k     =       \mathbf{r}_{k}^2 - d_{k}^2
          :label: eqnfepshakeconstr

with :math:`d_k` as before, so

.. math:: \begin{aligned}
          {\frac{\partial G}{\partial {\lambda}}}    &=&     -2 \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}
          :label: eqnfepshakeconstr2

Soft-core interactions: Beutler *et al.*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fig-softcore:

.. figure:: plots/softcore.*
   :height: 6.00000cm

   Soft-core interactions at :math:`{\lambda}=0.5`, with :math:`p=2` and
   :math:`C_6^A=C_{12}^A=C_6^B=C_{12}^B=1`.

In a free-energy calculation where particles grow out of nothing, or
particles disappear, using the simple linear interpolation of the
Lennard-Jones and Coulomb potentials as described in
:eq:`Equations %s <eqdVljdlambda>` and :eq:`%s <eqdVcoulombdlambda>` may lead to poor
convergence. When the particles have nearly disappeared, or are close to
appearing (at :math:`{\lambda}` close to 0 or 1), the interaction energy
will be weak enough for particles to get very close to each other,
leading to large fluctuations in the measured values of
:math:`\partial V/\partial {\lambda}` (which, because of the simple
linear interpolation, depends on the potentials at both the endpoints of
:math:`{\lambda}`).

To circumvent these problems, the singularities in the potentials need
to be removed. This can be done by modifying the regular Lennard-Jones
and Coulomb potentials with “soft-core” potentials that limit the
energies and forces involved at :math:`{\lambda}` values between 0 and
1, but not *at* :math:`{\lambda}=0` or 1.

In |Gromacs| the soft-core potentials :math:`V_{sc}` are shifted versions
of the regular potentials, so that the singularity in the potential and
its derivatives at :math:`r=0` is never reached. This formulation was
introduced by Beutler *et al.*\ :ref:`100 <refBeutler94>`:

.. math:: \begin{aligned}
          V_{sc}(r) &=& {(1-{\lambda})}V^A(r_A) + {\lambda}V^B(r_B)
              \\
          r_A &=& \left(\alpha \sigma_A^6 {\lambda}^p + r^6 \right)^\frac{1}{6}
              \\
          r_B &=& \left(\alpha \sigma_B^6 {(1-{\lambda})}^p + r^6 \right)^\frac{1}{6}\end{aligned}
          :label: eqnfepsoftcore

where :math:`V^A` and :math:`V^B` are the normal “hard core” Van der
Waals or electrostatic potentials in state A (:math:`{\lambda}=0`) and
state B (:math:`{\lambda}=1`) respectively, :math:`\alpha` is the
soft-core parameter (set with ``sc_alpha`` in the
:ref:`mdp` file), :math:`p` is the soft-core :math:`{\lambda}`
power (set with ``sc_power``), :math:`\sigma` is the radius
of the interaction, which is :math:`(C_{12}/C_6)^{1/6}` or an input
parameter (``sc_sigma``) when :math:`C_6` or :math:`C_{12}`
is zero. Beutler *et al.*\ :ref:`100 <refBeutler94>` probed various
combinations of the :math:`r` power values for the Lennard-Jones
and Coulombic interactions. |Gromacs| uses :math:`r^6` for both,
van der Waals and electrostatic interactions.

For intermediate :math:`{\lambda}`, :math:`r_A` and :math:`r_B` alter
the interactions very little for :math:`r > \alpha^{1/6} \sigma` and
quickly switch the soft-core interaction to an almost constant value for
smaller :math:`r` (:numref:`Fig. %s <fig-softcore>`). The force is:

.. math:: F_{sc}(r) = -\frac{\partial V_{sc}(r)}{\partial r} =
           {(1-{\lambda})}F^A(r_A) \left(\frac{r}{r_A}\right)^5 +
          {\lambda}F^B(r_B) \left(\frac{r}{r_B}\right)^5
          :label: eqnfepsoftcoreforce

where :math:`F^A` and :math:`F^B` are the “hard core” forces. The
contribution to the derivative of the free energy is:

.. math:: \begin{aligned}
          {\frac{\partial V_{sc}(r)}{\partial {\lambda}}} & = &
           V^B(r_B) -V^A(r_A)  + 
          	{(1-{\lambda})}\frac{\partial V^A(r_A)}{\partial r_A}
          		   \frac{\partial r_A}{\partial {\lambda}} + 
          	{\lambda}\frac{\partial V^B(r_B)}{\partial r_B}
          		   \frac{\partial r_B}{\partial {\lambda}}
          \nonumber\\
          &=&
           V^B(r_B) -V^A(r_A)  + \nonumber \\
           & &
           \frac{p \alpha}{6}
                 \left[ {\lambda}F^B(r_B) r^{-5}_B \sigma_B^6 {(1-{\lambda})}^{p-1} -
          	       {(1-{\lambda})}F^A(r_A) r^{-5}_A \sigma_A^6 {\lambda}^{p-1} \right]\end{aligned}
          :label: eqnfepsoftcorederivative

The original GROMOS Lennard-Jones soft-core
function\ :ref:`100 <refBeutler94>` uses :math:`p=2`, but :math:`p=1` gives a smoother
:math:`\partial H/\partial{\lambda}` curve. Another issue that should be
considered is the soft-core effect of hydrogens without Lennard-Jones
interaction. Their soft-core :math:`\sigma` is set with
``sc_sigma`` in the :ref:`mdp` file. These
hydrogens produce peaks in :math:`\partial H/\partial{\lambda}` at
:math:`{\lambda}` is 0 and/or 1 for :math:`p=1` and close to 0 and/or 1
with :math:`p=2`. Lowering ``sc_sigma``
will decrease this effect, but it will also increase the interactions
with hydrogens relative to the other interactions in the soft-core
state.

When soft-core potentials are selected (by setting
``sc_alpha >0``), and the Coulomb and Lennard-Jones
potentials are turned on or off sequentially, then the Coulombic
interaction is turned off linearly, rather than using soft-core
interactions, which should be less statistically noisy in most cases.
This behavior can be overwritten by using the :ref:`mdp` option
``sc-coul`` to ``yes``. Note that the
``sc-coul`` is only taken into account when lambda states
are used, not with ``couple-lambda0``  /
``couple-lambda1``, and you can still turn off soft-core
interactions by setting ``sc-alpha=0``. Additionally, the
soft-core interaction potential is only applied when either the A or B
state has zero interaction potential. If both A and B states have
nonzero interaction potential, default linear scaling described above is
used. When both Coulombic and Lennard-Jones interactions are turned off
simultaneously, a soft-core potential is used, and a hydrogen is being
introduced or deleted, the sigma is set to ``sc-sigma-min``,
which itself defaults to ``sc-sigma-default``.


Soft-core interactions: Gapsys *et al.*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we describe the functional form and parameters for 
the soft-cored non-bonded interactions using the formalism by Gapsys *et al.*\ :ref:`183 <refGapsys2012>`.

The Gapsys *et al.* soft-core is formulated to act on the level of van der Waals and electrostatic forces:
the non-bonded interactions are linearized at a point defined as, :math:`r_{scLJ}` or :math:`r_{scQ}`, respectively.
The linearization point depends on the state of the system as controlled by the :math:`\lambda` parameter and 
two parameters :math:`\alpha_Q` (set with ``sc-gapsys-scale-linpoint-q``) and :math:`\alpha_{LJ}` (set with ``sc-gapsys-scale-linpoint-lj``).
The dependence on :math:`\lambda` guarantees that the end-states are properly represented by their hard-core potentials.
:numref:`Fig. %s <fig-gapsyssc>` illustrates the behaviour of the linearization point, forces and integrated potential energies with respect
to the parameters :math:`\alpha_Q` and :math:`\alpha_{LJ}`. The optimal choices of the parameter values have been systematically explored in :ref:`183 <refGapsys2012>`. These recommended values are set by default when ``sc-function=gapsys`` is selected: ``sc-gapsys-scale-linpoint-q=0.3`` and ``sc-gapsys-scale-linpoint-lj=0.85``.

.. _fig-gapsyssc:

.. figure:: plots/gapsys-sc.*
        :width: 15.0cm

        Illustration of the soft-core parameter influence on the linearization point (top row), 
        forces (middle row) and energies (bottom row)
        for van der Waals (left column) and electrostatic interactions (right column).
        The case of two interacting atoms is considered.
        In state A both atoms have charges of 0.5 and :math:`\sigma=0.3` nm, :math:`\epsilon=0.5` kJ/mol.
        In state B all the non-bonded interactions are set to zero.
        The parameter :math:`\lambda` is set to 0.5 and electrostatic interaction cutoff is 1 nm.

The parameter :math:`\alpha_{LJ}` is a unitless scaling factor in the range :math:`[0,1)`.
It scales the position of the point from which the van der Waals force will be linearized.
The linearization of the force is allowed in the range :math:`[0,F_{min}^{LJ})`,
where setting :math:`\alpha_{LJ}=0` results in a standard hard-core van der Waals interaction.
Setting :math:`\alpha_{LJ}` closer to 1 brings the force linearization point towards 
the minimum in the Lennard-Jones force curve (:math:`F_{min}^{LJ}`).
This construct allows retaining the repulsion between two particles with non-zero C12 parameter at any :math:`\lambda` value.

The parameter :math:`\alpha_{Q}` has a unit of :math:`\frac{nm}{e^2}` and is defined in the range :math:`[0,\inf)`.
It scales the position of the point from which the Coulombic force will be linearized.
Even though in theory :math:`\alpha_{Q}` can be set to an arbitrarily large value,
algorithmically the linearization point for the force is bound in the range :math:`[0,F_{rcoul}^{Q})`,
where setting :math:`\alpha_{Q}=0` results in a standard hard-core Coulombic interaction.
Setting :math:`\alpha_{Q}` to a larger value softens the Coulombic force.

In all the notations below, for simplicity, the distance between two atoms :math:`i` and :math:`j` is noted as :math:`r`, i.e. :math:`r=r_{ij}`.

Forces: van der Waals interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: \begin{aligned}
          \mathbf{F}_{ij}^{LJ}(\mathbf{r})=\begin{cases}
          (\frac{12C_{ij}^{(12)}}{r^{13}} - \frac{6C_{ij}^{(6)}}{r^7})\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scLJ}$} 
          \\
          \frac{d\mathbf{F}_{ij}^{LJ}}{dr}_{r=r_{scLJ}}r + \mathbf{F}_{ij}^{LJ}(r_{scLJ}), & \mbox{if } \mbox{ $r<r_{scLJ}$}
          \end{cases}\end{aligned}
          :label: eqvdwforces

where the switching point between the soft and hard-core Lennard-Jones forces
:math:`r_{scLJ} = \alpha_{LJ}(\frac{26}{7}\sigma^6\lambda)^{\frac{1}{6}}` for state A, and
:math:`r_{scLJ} = \alpha_{LJ}(\frac{26}{7}\sigma^6(1-\lambda))^{\frac{1}{6}}` for state B.
In analogy to the Beutler *et al.* soft core version, :math:`\sigma` is the radius of the interaction, which is :math:`(C_{12}/C_6)^{1/6}`
or an input parameter (set with ``sc-sigma-LJ-gapsys``) when C6 or C12 is zero. The default value for this parameter is ``sc-sigma-LJ-gapsys=0.3``.

Explicit expression:

.. math:: \begin{aligned}
          \mathbf{F}_{LJ}(\mathbf{r})=\begin{cases}
          \left(\frac{12C^{(12)}}{r^{13}} - \frac{6C^{(6)}}{r^7}\right)\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scLJ}$} 
          \\
          \left(-\frac{156C^{(12)}}{r_{scLJ}^{14}} + \frac{42C^{(6)}}{r_{scLJ}^{8}}\right)\mathbf{r} + \frac{168C^{(12)}}{r_{scLJ}^{13}} - \frac{48C^{(6)}}{r_{scLJ}^{7}}, & \mbox{if } \mbox{ $r<r_{scLJ}$}
          \end{cases}\end{aligned}
          :label: eqvdwforcesexpl

Forces: Coulomb interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: \begin{aligned}
          \mathbf{F}_{ij}^{Q}(\mathbf{r})=\begin{cases}
          \frac{q_{i}q_{j}}{4{\pi}{\varepsilon_0}{\varepsilon_r}r^{2}}\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scQ} < r_{cutoffQ}$} 
          \\
          \frac{d\mathbf{F}_{ij}^{Q}}{dr}_{r=r_{scQ}}r + \mathbf{F}_{ij}^{Q}(r_{scQ}), & \mbox{if } \mbox{ $r<r_{scQ} < r_{cutoffQ}$}
          \\
          \frac{d\mathbf{F}_{ij}^{Q}}{dr}_{r=r_{cutoffQ}}r + \mathbf{F}_{ij}^{Q}(r_{cutoffQ}), & \mbox{if } \mbox{ $r < r^{scQ} \geq r_{cutoffQ}$} 
          \end{cases}\end{aligned}
          :label: eqqforces

where the switching point :math:`r^{sc}` between the soft and hard-core electrostatic forces is 
:math:`r_{scQ} = \alpha_Q(1+|q_iq_j|)\lambda^{\frac{1}{6}}` for state A, and
:math:`r_{scQ} = \alpha_Q(1+|q_iq_j|)(1-\lambda)^{\frac{1}{6}}` for state B.
The :math:`\lambda` dependence of the linearization point for both van der Waals and Coulombic interactions is of the same power :math:`1/6`.

Explicit expression:

.. math:: \begin{aligned}
          \mathbf{F}_{Q}(\mathbf{r})=\begin{cases}
          \frac{q_iq_j}{4{\pi}{\varepsilon_0}{\varepsilon_r}r^{2}}\frac{\mathbf{r}}{r}, & \mbox{if } \mbox{ $r \geq r_{scQ} < r_{cutoffQ}$} 
          \\
          \frac{1}{4{\pi}{\varepsilon_0}{\varepsilon_r}}\big( -\frac{2q_{i}q_{j}}{r_{sc}^3}\mathbf{r} + \frac{3q_iq_j}{r_{sc}^2} \big), & \mbox{if } \mbox{ $r<r_{scQ} < r_{cutoffQ}$}
          \\
          \frac{1}{4{\pi}{\varepsilon_0}{\varepsilon_r}}\big( -\frac{2q_{i}q_{j}}{r_{cutoffQ}^3}\mathbf{r} + \frac{3q_iq_j}{r_{cutoffQ}^2} \big), & \mbox{if } \mbox{ $r < r_{scQ} \geq r_{cutoffQ}$}                        \end{cases}\end{aligned}
          :label: eqqforcesexpl

Energies: van der Waals interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Explicition definition of energies:

.. math:: \begin{aligned}
          V_{LJ}(r)=\begin{cases}
          \frac{C^{(12)}}{r^{12}} - \frac{C^{(6)}}{r^6}, & \mbox{if } \mbox{ $r \geq r_{scLJ}$} 
          \\
          \left(\frac{78C^{(12)}}{r_{scLJ}^{14}} - \frac{21C^{(6)}}{r_{scLJ}^{8}}\right)r^2 - \left(\frac{168C^{(12)}}{r_{scLJ}^{13}} - \frac{48C^{(12)}}{r_{scLJ}^{7}}\right)r
          + \frac{91C^{(12)}}{r_{scLJ}^{12}} - \frac{28C^{(6)}}{r_{scLJ}^{6}}, & \mbox{if } \mbox{ $r<r_{scLJ}$}
          \end{cases}\end{aligned}
          :label: eqvdwener

Energies: Coulomb interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: \begin{aligned}
          V_{Q}(r)=\begin{cases}
          \frac{q_{i}q_{j}}{4{\pi}{\varepsilon_0}{\varepsilon_r}r}, & \mbox{if } \mbox{ $r \geq r_{scQ} < r_{cutoffQ}$}
          \\
          \frac{q_{i}q_{j}}{r_{scQ}^3}r^2 - \frac{3q_iq_j}{r_{scQ}^2}r + \frac{3q_iq_j}{r_{scQ}}, & \mbox{if } \mbox{ $r<r_{scQ} < r_{cutoffQ}$}
          \\
          \frac{q_{i}q_{j}}{r_{cutoffQ}^3}r^2 - \frac{3q_iq_j}{r_{cutoffQ}^2}r + \frac{3q_iq_j}{r_{cutoffQ}}, & \mbox{if } \mbox{ $r < r_{scQ} \geq r_{cutoffQ}$}                
          \end{cases}\end{aligned}
          :label: eqqener

:math:`\partial H / \partial \lambda`: van der Waals interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we provide the explicit expressions of :math:`\partial H/ \partial \lambda` for Lennard-Jones potential, when :math:`r<r_{scLJ}`.
For simplicity, in the expression below we use the notation :math:`r_{scLJ_A}=r_{scA}` and :math:`r_{scLJ_B}=r_{scB}`.

.. math:: \begin{aligned}
          \frac{\partial{H}}{\partial{\lambda}} &= V_{LJ}^B(r) - V_{LJ}^A(r) + (1-\lambda)\frac{\partial{V_{LJ}^A(r)}}{\partial{\lambda}} + \lambda\frac{\partial{V_{LJ}^B(r)}}{\partial{\lambda}} \\
          & =  \left(\frac{78C^{(12)}_B}{r_{scB}^{14}} - \frac{21C^{(6)}_B}{r_{scB}^{8}}\right)r^2 - \left(\frac{168C^{(12)}_B}{r_{scB}^{13}} - \frac{48C^{(12)}_B}{r_{scB}^{7}}\right)r
          + \frac{91C^{(12)}_B}{r_{scB}^{12}} - \frac{28C^{(6)}_B}{r_{scB}^{6}} \\
          & -  \left[\left(\frac{78C^{(12)}_A}{r_{scA}^{14}} - \frac{21C^{(6)}_A}{r_{scA}^{8}}\right)r^2 - \left(\frac{168C^{(12)}_A}{r_{scA}^{13}} - \frac{48C^{(12)}_A}{r_{scA}^{7}}\right)r
          + \frac{91C^{(12)}_A}{r_{scA}^{12}} - \frac{28C^{(6)}_A}{r_{scA}^{6}} \right]\\
          & +  \frac{14(\lambda-1)}{\lambda}\left[\left(\frac{13C^{(12)}_A}{r_{scA}^{14}} - \frac{2C^{(6)}_A}{r_{scA}^{8}}\right)r^2
          - \left(\frac{26C^{(12)}_A}{r_{scA}^{13}} - \frac{4C^{(6)}_A}{r_{scA}^{7}}\right)r
          + \frac{13C^{(12)}_A}{r_{scA}^{12}} - \frac{2C^{(6)}_A}{r_{scA}^{6}}\right] \\
          & +  \frac{14\lambda}{1-\lambda}\left[\left(\frac{13C^{(12)}_B}{r_{scB}^{14}} - \frac{2C^{(6)}_B}{r_{scB}^{8}}\right)r^2
          - \left(\frac{26C^{(12)}_B}{r_{scB}^{13}} - \frac{4C^{(6)}_B}{r_{scB}^{7}}\right)r
          + \frac{13C^{(12)}_B}{r_{scB}^{12}} - \frac{2C^{(6)}_B}{r_{scB}^{6}}\right] \end{aligned}
          :label: eqvdwdhdl

:math:`\partial H/ \partial \lambda` for Lennard-Jones potential, when :math:`r \geq r_{scLJ}` is calculated as a standard hard-core
contribution to :math:`\partial H/ \partial \lambda`: :math:`\frac{\partial{H}}{\partial{\lambda}} = V_{LJ}^B(r) - V_{LJ}^A(r)`.

:math:`\partial H/ \partial \lambda` for Coulomb interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we provide the explicit expressions of :math:`\partial H/ \partial \lambda` for Coulomb potential, when :math:`r<r_{scQ}<r_{cutoffQ}`.
For simplicity, in the expression below we use the notation :math:`r_{scQ_A}=r_{scA}` and :math:`r_{scQ_B}=r_{scB}`.

.. math:: \begin{aligned}
          \frac{\partial{H}}{\partial{\lambda}} &= V_Q^B(r) - V_Q^A(r) + (1-\lambda)\frac{\partial{V_Q^A(r)}}{\partial{\lambda}} + \lambda\frac{\partial{V_Q^B(r)}}{\partial{\lambda}} \\
          & =  \frac{q_{i}^Bq_{j}^B}{r_{scB}^3}r^2 - \frac{3q_i^Bq_j^B}{r_{scB}^2}r + \frac{3q^B_iq_j^B}{r_{scB}} \\
          & -  \left[\frac{q_{i}^Aq_{j}^A}{r_{scA}^3}r^2 - \frac{3q_i^Aq_j^A}{r_{scA}^2}r + \frac{3q^A_iq_j^A}{r_{scA}}\right] \\
          & +  \frac{\lambda-1}{2\lambda}\left[\frac{q_i^Aq_j^A}{r_{scA}^3}r^2 - \frac{2q_i^Aq_j^A}{r_{scA}^2}r + \frac{q_i^Aq_j^A}{r_{scA}}\right] \\
          & +  \frac{\lambda}{2(1-\lambda)}\left[\frac{q_i^Bq_j^B}{r_{scB}^3}r^2 - \frac{2q_i^Bq_j^B}{r_{scB}^2}r + \frac{q_i^Bq_j^B}{r_{scB}}\right] \end{aligned}
          :label: eqqdhdl

:math:`\partial H/ \partial \lambda` for Coulomb potential, when :math:`r < r_{scQ} \geq r_{cutoffQ}` is calculated using the same expression above
by setting :math:`r_{scA}=r_{cutoffQ}` and :math:`r_{scB}=r_{cutoffQ}`.

:math:`\partial H/ \partial \lambda` for Coulomb potential, when :math:`r \geq r_{scQ} < r_{cutoffQ}` is calculated as a standard hard-core
contribution to :math:`\partial H/ \partial \lambda`: :math:`\frac{\partial{H}}{\partial{\lambda}} = V_{Q}^B(r) - V_{Q}^A(r)`.



