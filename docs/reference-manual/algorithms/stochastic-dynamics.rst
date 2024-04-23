Stochastic Dynamics
-------------------

Stochastic or velocity Langevin dynamics adds a friction and a noise
term to Newton’s equations of motion, as

.. math::  m_i {{\mbox{d}}^2 \mathbf{r}_i \over {\mbox{d}}t^2} =
   - m_i \gamma_i {{\mbox{d}}\mathbf{r}_i \over {\mbox{d}}t} + \mathbf{F}_i(\mathbf{r}) + {\stackrel{\circ}{\mathbf{r}}}_i,
   :label: eqnSDeq

where :math:`\gamma_i` is the friction constant :math:`[1/\mbox{ps}]`
and :math:`{\stackrel{\circ}{\mathbf{r}}}_i\!\!(t)` is a
noise process with
:math:`\langle {\stackrel{\circ}{r}}_i\!\!(t) {\stackrel{\circ}{r}}_j\!\!(t+s) \rangle = 2 m_i \gamma_i k_B T \delta(s) \delta_{ij}`. When :math:`1/\gamma_i`
is large compared to the time scales present in the system, one could
see stochastic dynamics as molecular dynamics with stochastic
temperature-coupling. But any processes that take longer than
:math:`1/\gamma_i`, e.g. hydrodynamics, will be dampened. Since each
degree of freedom is coupled independently to a heat bath, equilibration
of fast modes occurs rapidly. For simulating a system in vacuum there is
the additional advantage that there is no accumulation of errors for the
overall translational and rotational degrees of freedom. When
:math:`1/\gamma_i` is small compared to the time scales present in the
system, the dynamics will be completely different from MD, but the
sampling is still correct.

In |Gromacs| there is one simple and efficient implementation. Its
accuracy is equivalent to the normal MD leap-frog and Velocity Verlet
integrator. It is nearly identical to the common way of discretizing the
Langevin equation, but the friction and velocity term are applied in an
impulse fashion \ :ref:`51 <refGoga2012>`. It can be described as:

.. math::  \begin{aligned}
   \mathbf{v}'  &~=~&   \mathbf{v}(t-{{\frac{1}{2}}{{\Delta t}}}) + \frac{1}{m}\mathbf{F}(t){{\Delta t}}\\
   \Delta\mathbf{v}     &~=~&   -\alpha \, \mathbf{v}'(t+{{\frac{1}{2}}{{\Delta t}}}) + \sqrt{\frac{k_B T}{m} \alpha (2 - \alpha)} \, {\mathbf{r}^G}_i \\
   \mathbf{r}(t+{{\Delta t}})   &~=~&   \mathbf{r}(t)+\left(\mathbf{v}' +\frac{1}{2}\Delta \mathbf{v}\right){{\Delta t}}
  \end{aligned}
  :label: eqnsd1int

.. math:: \begin{aligned}
   \mathbf{v}(t+{{\frac{1}{2}}{{\Delta t}}})  &~=~&   \mathbf{v}' + \Delta \mathbf{v} \\
   \alpha &~=~& 1 - e^{-\gamma {{\Delta t}}}\end{aligned}
   :label: eqnsd1xupd

where :math:`{\mathbf{r}^G}_i` is Gaussian distributed
noise with :math:`\mu = 0`, :math:`\sigma = 1`. The velocity is first
updated a full time step without friction and noise to get
:math:`\mathbf{v}'`, identical to the normal update in
leap-frog. The friction and noise are then applied as an impulse at step
:math:`t+{{\Delta t}}`. The advantage of this scheme is that the
velocity-dependent terms act at the full time step, which makes the
correct integration of forces that depend on both coordinates and
velocities, such as constraints and dissipative particle dynamics (DPD,
not implemented yet), straightforward. With constraints, the coordinate
update :eq:`eqn. %s <eqnsd1xupd>` is split into a normal leap-frog update
and a :math:`\Delta \mathbf{v}`. After both of these
updates the constraints are applied to coordinates and velocities.

When using SD as a thermostat, an appropriate value for :math:`\gamma`
is e.g. 0.5 ps\ :math:`^{-1}`, since this results in a friction that is
lower than the internal friction of water, while it still provides
efficient thermostatting.
