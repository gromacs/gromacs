.. _awh:

Adaptive biasing with AWH
-------------------------

The accelerated weight histogram method :ref:`185 <refLidmar2012>`
:ref:`137 <reflindahl2014accelerated>` calculates the PMF along a reaction coordinate by adding
an adaptively determined biasing potential. AWH flattens free energy
barriers along the reaction coordinate by applying a history-dependent
potential to the system that “fills up” free energy minima. This is
similar in spirit to other adaptive biasing potential methods, e.g. the
Wang-Landau \ :ref:`138 <refwang2001efficient>`, local
elevation \ :ref:`139 <refhuber1994local>` and
metadynamics \ :ref:`140 <reflaio2002escaping>` methods.
The initial sampling stage of AWH makes the method robust against the
choice of input parameters. Furthermore, the target distribution along
the reaction coordinate may be chosen freely.

Basics of the method
^^^^^^^^^^^^^^^^^^^^

The AWH method can act on two different types of reaction coordinates.
It can work directly on a discrete reaction coordinate :math:`\lambda`
in case this is the free-energy coupling parameter :ref:`187 <reflundborg2021>`.
And it can act on reaction coordinates that are (continuous) functions of the coordinates:
:math:`\xi(x)`. In this case AWH acts on a *reference coordinate* :math:`\lambda` which
takes discrete values and is coupled to :math:`\xi(x)` using an umbrella function :math:`Q`.
We will now describe the method for the most general case. When acting directly
on :math:`\lambda`, the function :math:`Q` is zero. The fundamentals of the
method are based on the connection between atom coordinates and :math:`\lambda` and
are established through the extended ensemble :ref:`68 <refLyubartsev1992>`,

.. math:: P(x,\lambda) = \frac{1}{\mathcal{Z}}e^{g(\lambda) - Q(\xi(x),\lambda) - V(x)},
          :label: eqawhpxlambda

where :math:`g(\lambda)` is a bias function (a free variable) and
:math:`V(x)` is the unbiased potential energy of the system. The
distribution along :math:`\lambda` can be tuned to be any predefined
*target distribution* :math:`\rho(\lambda)` (often chosen to be flat) by
choosing :math:`g(\lambda)` wisely. This is evident from

.. math:: P(\lambda) = \int P(x,\lambda)  dx =
          \frac{1}{\mathcal{Z}}e^{g(\lambda)} \int e^{- Q(\xi(x),\lambda) - V(x)}  dx
          \equiv \frac{1}{\mathcal{Z}}e^{g(\lambda) - F(\lambda)},
          :label: eqawhplambda

where :math:`F(\lambda)` is the free energy

.. math:: F(\lambda) = -\ln \int e^{- Q(\xi(x),\lambda) - V(x)}  dx.
          :label: eqawhflambda

The reaction coordinate :math:`\xi(x)` is commonly coupled to
:math:`\lambda` with a harmonic potential

.. math:: Q(\xi,\lambda) = \frac{1}{2} \beta k (\xi - \lambda)^2,
          :label: eqnawhbasic

so that for large force constants :math:`k`,
:math:`\xi \approx \lambda`. Note the use of dimensionless energies for
compatibility with previously published work. Units of energy are
obtained by multiplication with :math:`k_BT=1/\beta`. In the simulation,
:math:`\lambda` samples the user-defined sampling interval :math:`I`.

Being the convolution of the PMF with the Gaussian defined by the
harmonic potential, :math:`F(\lambda)` is a smoothened version of the
PMF. :eq:`Eq. %s <eqawhplambda>` shows that in order to obtain
:math:`P(\lambda)=\rho(\lambda)`, :math:`F(\lambda)` needs to be
determined accurately. Thus, AWH adaptively calculates
:math:`F(\lambda)` and simultaneously converges :math:`P(\lambda)`
toward :math:`\rho(\lambda)`.

For a multidimensional reaction coordinate :math:`\xi`, the sampling
interval is the Cartesian product :math:`I=\Pi_{\mu} I_{\mu}` (a rectangular
domain).

N.b., it is not yet possible to use AWH for alchemical transformations
that involve perturbed masses or constraints.

The free energy update
^^^^^^^^^^^^^^^^^^^^^^

AWH is initialized with an estimate of the free energy
:math:`F_0(\lambda)`. At regular time intervals this estimate is updated
using data collected in between the updates. At update :math:`n`, the
applied bias :math:`g_n(\lambda)` is a function of the current free
energy estimate :math:`F_n(\lambda)` and target distribution
:math:`\rho_n(\lambda)`,

.. math:: g_n(\lambda) = \ln \rho_n(\lambda) +F_n(\lambda),
          :label: eqawhgrhofrelation

which is consistent with :eq:`Eq. %s <eqawhplambda>`. Note that also the
target distribution may be updated during the simulation (see examples
in section :ref:`awhtargets`). Substituting this choice of :math:`g=g_n`
back into :eq:`Eq. %s <eqawhplambda>` yields the simple free energy update

.. math:: \Delta F_n(\lambda)
          = F(\lambda) - F_n(\lambda)
          = -\ln\frac{P_n(\lambda)}{\rho_n(\lambda)},
          :label: eqawhdfnaive

which would yield a better estimate :math:`F_{n+1} = F_n + \Delta F_n`,
assuming :math:`P_n(\lambda)` can be measured accurately. AWH estimates
:math:`P_n(\lambda)` by regularly calculating the conditional
distribution

.. math:: \omega_n(\lambda|x) \equiv P_n(\lambda|x) = \frac{e^{g_n(\lambda) - Q(\xi(x), \lambda)}}{\sum_{\lambda'} e^{g_n(\lambda') - Q(\xi(x),\lambda')}}.
          :label: eqawhomega

Accumulating these probability weights yields
:math:`\sum_t \omega(\lambda|x(t)) \sim P_n(\lambda)`, where
:math:`\int P_n(\lambda|x) P_n(x) dx = P_n(\lambda)` has been used. The
:math:`\omega_n(\lambda|x)` weights are thus the samples of the AWH
method. With the limited amount of sampling one has in practice, update
scheme :eq:`%s <eqawhdfnaive>` yields very noisy results. AWH instead applies a
free energy update that has the same form but which can be applied
repeatedly with limited and localized sampling,

.. math:: \Delta F_n = -\ln \frac{W_n(\lambda) + \sum_t \omega_n(\lambda|x(t))}{W_n(\lambda) + \sum_t\rho_n(\lambda)) }.
          :label: eqnawhsampling

Here :math:`W_n(\lambda)` is the *reference weight histogram*
representing prior sampling. The update for :math:`W(\lambda)`,
disregarding the initial stage (see section :ref:`awhinitialstage`), is

.. math:: W_{n+1}(\lambda) = W_n(\lambda) + \sum_t\rho_n(\lambda).
          :label: eqawhwupdate

Thus, the weight histogram equals the targeted, “ideal” history of
samples. There are two important things to note about the free energy
update. First, sampling is driven away from oversampled, currently local
regions. For such :math:`\lambda` values,
:math:`\omega_n(\lambda) > \rho_n(\lambda)` and
:math:`\Delta F_n(\lambda) < 0`, which by :eq:`Eq. %s <eqawhgrhofrelation>`
implies :math:`\Delta g_n(\lambda) < 0` (assuming
:math:`\Delta \rho_n \equiv 0`). Thus, the probability to sample
:math:`\lambda` decreases after the update (see :eq:`Eq. %s <eqawhplambda>`).
Secondly, the normalization of the histogram
:math:`N_n=\sum_\lambda W_n(\lambda)`, determines the update size
:math:`| \Delta F(\lambda) |`. For instance, for a single sample
:math:`\omega(\lambda|x)`, and using a harmonic potential
(:see :eq:`Eq. %s <eqnawhbasic>`),
the shape of the update is approximately a
Gaussian function of width :math:`\sigma=1/\sqrt{\beta k}` and height
:math:`\propto 1/N_n` :ref:`137 <reflindahl2014accelerated>`,

.. math:: | \Delta F_n(\lambda) | \propto \frac{1}{N_n} e^{-\frac{1}{2} \beta k (\xi(x) - \lambda)^2}.
          :label: eqawhdfsize

When directly controlling the lambda state of the system, the shape of
the update is instead

.. math:: | \Delta F_n(\lambda) | \propto \frac{1}{N_n} P_n(\lambda | x).
          :label: eqawhdfsizelambda

Therefore, in both cases, as samples accumulate in :math:`W(\lambda)` and
:math:`N_n` grows, the updates get smaller, allowing for the free energy to
converge.

Note that quantity of interest to the user is not :math:`F(\lambda)` but
the PMF :math:`\Phi(\xi)`. :math:`\Phi(\xi)` is extracted by reweighting
samples :math:`\xi(t)` on the fly \ :ref:`137 <reflindahl2014accelerated>` (see
also section :ref:`awhreweight`) and will converge at the same rate as
:math:`F(\lambda)`, see :numref:`Fig. %s <fig-awhbiasevolution1>`. The PMF will be
written to output (see section :ref:`awhusage`).

Applying the bias to the system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bias potential can be applied to the system in two ways. Either by
applying a harmonic potential centered at :math:`\lambda(t)`, which is
sampled using (rejection-free) Monte-Carlo sampling from the conditional
distribution :math:`\omega_n(\lambda | x(t)) = P_n(\lambda | x(t))`, see
:eq:`Eq. %s <eqawhomega>`. This is also called Gibbs sampling or independence
sampling. Alternatively, and by default in the code, the following
*convolved bias potential* can be applied,

.. math:: U_n(\xi) = -\ln \int e^{ g_n(\lambda) -Q(\xi,\lambda)} d \lambda.
          :label: eqawhbiaspotential

These two approaches are equivalent in the sense that they give rise to
the same biased probabilities :math:`P_n(x)`
(cf. :eq:`%s <eqawhpxlambda>`) while the dynamics are clearly
different in the two cases. This choice does not affect the internals of
the AWH algorithm, only what force and potential AWH returns to the MD
engine.

Along a bias dimension directly controlling the
:math:`\lambda` state of the system, such as when controlling free energy
perturbations, the Monte-Carlo sampling alternative is always used, even if
a convolved bias potential is chosen to be used along the other dimensions
(if there are more than one).


.. _fig-awhbiasevolution1:

.. figure:: plots/awh-traj.*
        :width: 8.0cm

        AWH evolution in time for a Brownian particle in a double-well
        potential. The reaction coordinate :math:`\xi(t)` traverses the sampling
        interval multiple times in the initial stage before exiting and entering
        the final stage. In the final stage, the dynamics of
        :math:`\xi` becomes increasingly diffusive.

.. _fig-awhbiasevolution2:

.. figure:: plots/awh-invN.*
        :width: 8.0cm

        In the final stage, the dynamics of
        :math:`\xi` becomes increasingly diffusive. The times of covering are
        shown as :math:`\times`-markers of different colors. At these times the
        free energy update size :math:`\sim 1/N`, where :math:`N` is the size of
        the weight histogram, is decreased by scaling :math:`N` by a factor of
        :math:`\gamma=3` (note that the default value of :math:`\gamma`
        is 2 since |Gromacs| 2024).

.. _fig-awhbiasevolution3:

.. figure:: plots/awh-sampleweights.*
        :width: 8.0cm

        In the final stage, :math:`N` grows at the
        sampling rate and thus :math:`1/N\sim1/t`. The exit from the final stage
        is determined on the fly by ensuring that the effective sample weight
        :math:`s` of data collected in the final stage exceeds that of initial
        stage data (note that :math:`\ln s(t)` is plotted).

.. _fig-awhbiasevolution4:

.. figure:: plots/awh-pmfs.*
        :width: 8.0cm

        An estimate of the PMF is also extracted from the simulation (bottom
        right), which after exiting the initial stage should estimate global
        free energy differences fairly accurately.

.. _awhinitialstage:

The initial stage
~~~~~~~~~~~~~~~~~

Initially, when the bias potential is far from optimal, samples will be
highly correlated. In such cases, letting :math:`W(\lambda)` accumulate
samples as prescribed by :eq:`Eq. %s <eqawhwupdate>`, entails
a too rapid decay of the free energy update size. This motivates
splitting the simulation into an *initial stage* where the weight
histogram grows according to a more restrictive and robust protocol, and
a *final stage* where the weight histogram grows linearly at the
sampling rate (:eq:`Eq. %s <eqawhwupdate>`). The AWH initial
stage takes inspiration from the well-known Wang-Landau algorithm \ :ref:`138 <refwang2001efficient>`,
although there are differences in the details.

In the initial stage the update size is kept constant (by keeping
:math:`N_n` constant) until a transition across the sampling interval
has been detected, a “covering”. For the definition of a covering, see
:eq:`Eq. %s <eqawhcovering>` below. After a covering has
occurred, :math:`N_n` is scaled up by a constant “growth factor”
:math:`\gamma`, which by default has the value of 2. Thus, in the
initial stage :math:`N_n` is set dynamically as
:math:`N_{n} = \gamma^{m} N_0`, where :math:`m` is the number of
coverings. Since the update size scales as :math:`1/N` (
:eq:`Eq. %s <eqawhdfsize>`) this leads to a close to
exponential decay of the update size in the initial stage, see
:numref:`Fig. %s <fig-awhbiasevolution1>`.

The update size directly determines the rate of change of
:math:`F_n(\lambda)` and hence, from
:eq:`Eq. %s <eqawhgrhofrelation>`, also the rate of change of
the bias funcion :math:`g_n(\lambda)` Thus initially, when :math:`N_n`
is kept small and updates large, the system will be driven along the
reaction coordinate by the constantly fluctuating bias. If :math:`N_0`
is set small enough, the first transition will typically be fast because
of the large update size and will quickly give a first rough estimate of
the free energy. The second transition, using :math:`N_1=\gamma N_0`
refines this estimate further. Thus, rather than very carefully filling
free energy minima using a small initial update size, the sampling
interval is sweeped back-and-forth multiple times, using a wide range of
update sizes, see :numref:`Fig. %s <fig-awhbiasevolution1>`. This
way, the initial stage also makes AWH robust against the choice of
:math:`N_0`.

The covering criterion
^^^^^^^^^^^^^^^^^^^^^^

In the general case of a multidimensional reaction coordinate
:math:`\lambda=(\lambda_{\mu})`, the sampling interval :math:`I` is
considered covered when all dimensions have been covered. A dimension
:math:`d` is covered if all points :math:`\lambda_{\mu}` in the
one-dimensional sampling interval :math:`I_{\mu}` have been “visited”.
Finally, a point :math:`\lambda_{\mu} \in I_{\mu}` has been visited if there is
at least one point :math:`\lambda^*\in I` with
:math:`\lambda^*_{\mu} = \lambda_{\mu}` that since the last covering has
accumulated probability weight corresponding to the peak of a
multidimensional Gaussian distribution

.. math:: \Delta W(\lambda^*)
          \ge w_{\mathrm{peak}}
          \equiv \prod_{\mu} \frac{\Delta \lambda_{\mu}}{\sqrt{2\pi}\sigma_k}.
          :label: eqawhcovering

Here, :math:`\Delta \lambda_{\mu}` is the point spacing of the discretized
:math:`I_{\mu}` and :math:`\sigma_k=1/\sqrt{\beta k_{\mu}}` (where :math:`k_{\mu}`
is the force constant) is the Gaussian width.

Exit from the initial stage
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For longer times, when major free energy barriers have largely been
flattened by the converging bias potential, the histogram
:math:`W(\lambda)` should grow at the actual sampling rate and the
initial stage needs to be exited \ :ref:`141 <refbelardinelli2007fast>`.
There are multiple reasonable (heuristic) ways of determining when this
transition should take place. One option is to postulate that the number
of samples in the weight histogram :math:`N_n` should never exceed the
actual number of collected samples, and exit the initial stage when this
condition breaks \ :ref:`137 <reflindahl2014accelerated>`. In the initial stage,
:math:`N` grows close to exponentially while the collected number of
samples grows linearly, so an exit will surely occur eventually. Here we
instead apply an exit criterion based on the observation that
“artificially” keeping :math:`N` constant while continuing to collect
samples corresponds to scaling down the relative weight of old samples
relative to new ones. Similarly, the subsequent scaling up of :math:`N`
by a factor :math:`\gamma` corresponds to scaling up the weight of old
data. Briefly, the exit criterion is devised such that the weight of a
sample collected *after* the initial stage is always larger or equal to
the weight of a sample collected *during* the initial stage, see
:numref:`Fig. %s <fig-awhbiasevolution1>`. This is consistent with
scaling down early, noisy data.

The initial stage exit criterion will now be described in detail. We
start out at the beginning of a covering stage, so that :math:`N` has
just been scaled by :math:`\gamma` and is now kept constant. Thus, the
first sample of this stage has the weight :math:`s= 1/\gamma` relative
to the last sample of the previous covering stage. We assume that
:math:`\Delta N` samples are collected and added to :math:`W` for each
update . To keep :math:`N` constant, :math:`W` needs to be scaled down
by a factor :math:`N/(N + \Delta N)` after every update. Equivalently,
this means that new data is scaled up relative to old data by the
inverse factor. Thus, after :math:`\Delta n` updates a new sample has
the relative weight
:math:`s=(1/\gamma) [(N_n + \Delta N)/N_n]^{\Delta n}`. Now assume
covering occurs at this time. To continue to the next covering stage,
:math:`N` should be scaled by :math:`\gamma`, which corresponds to again
multiplying :math:`s` by :math:`1/\gamma`. If at this point
:math:`s \ge \gamma`, then after rescaling :math:`s \ge 1`; i.e. overall
the relative weight of a new sample relative to an old sample is still
growing fast. If on the contrary :math:`s < \gamma`, and this defines
the exit from the initial stage, then the initial stage is over and from
now :math:`N` simply grows at the sampling rate (see
:eq:`Eq. %s <eqawhwupdate>`). To really ensure that
:math:`s\ge 1` holds before exiting, so that samples after the exit have
at least the sample weight of older samples, the last covering stage is
extended by a sufficient number of updates.

.. _awhtargets:

Choice of target distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The target distribution :math:`\rho(\lambda)` is traditionally chosen to
be uniform

.. math:: \rho_{\mathrm{const}}(\lambda) = \mathrm{const.}
          :label: eqnawhuniformdist

This choice exactly flattens :math:`F(\lambda)` in user-defined
sampling interval :math:`I`. Generally,
:math:`\rho(\lambda)=0, \lambda\notin I`. In certain cases other choices
may be preferable. For instance, in the multidimensional case the
rectangular sampling interval is likely to contain regions of very high
free energy, e.g. where atoms are clashing. To exclude such regions,
:math:`\rho(\lambda)` can specified by the following function of the
free energy

.. math:: \rho_{\mathrm{cut}}(\lambda) \propto \frac{1}{1+ e^{F(\lambda) - F_{\mathrm{cut}}}},
          :label: eqawhrhocut


where :math:`F_{\mathrm{cut}}` is a free energy cutoff (relative to
:math:`\min_\lambda F(\lambda)`). Thus, regions of the sampling interval
where :math:`F(\lambda) > F_{\mathrm{cut}}` will be exponentially
suppressed (in a smooth fashion). Alternatively, very high free energy
regions could be avoided while still flattening more moderate free
energy barriers by targeting a Boltzmann distribution corresponding to
scaling :math:`\beta=1/k_BT` by a factor :math:`0<s_\beta<1`,

.. math:: \rho_{\mathrm{Boltz}}(\lambda) \propto e^{-s_\beta F(\lambda)},
          :label: eqawhrhoboltz

The parameter :math:`s_\beta` determines to what degree the free energy
landscape is flattened; the lower :math:`s_\beta`, the flatter. Note
that both :math:`\rho_{\mathrm{cut}}(\lambda)` and
:math:`\rho_{\mathrm{Boltz}}(\lambda)` depend on :math:`F(\lambda)`,
which needs to be substituted by the current best estimate
:math:`F_n(\lambda)`. Thus, the target distribution is also updated
(consistently with :eq:`Eq. %s <eqawhgrhofrelation>`).

There is in fact an alternative approach to obtaining
:math:`\rho_{\mathrm{Boltz}}(\lambda)` as the limiting target
distribution in AWH, which is particular in the way the weight histogram
:math:`W(\lambda)` and the target distribution :math:`\rho` are updated
and coupled to each other. This yields an evolution of the bias
potential which is very similar to that of well-tempered
metadynamics \ :ref:`142 <refbarducci2008well>`,
see \ :ref:`137 <reflindahl2014accelerated>` for details. Because of the popularity and
success of well-tempered metadynamics, this is a special case worth
considering. In this case :math:`\rho` is a function of the reference
weight histogram

.. math:: \rho_{\mathrm{Boltz,loc}}(\lambda) \propto W(\lambda),
          :label: eqnawhweighthistogram

and the update of the weight histogram is modified (cf.
:eq:`Eq. %s <eqawhwupdate>`)

.. math:: W_{n+1}(\lambda) =  W_{n}(\lambda) + s_{\beta}\sum_t \omega(\lambda | x(t)).
          :label: eqnawhupdateweighthist

Thus, here the weight histogram equals the real history of samples, but
scaled by :math:`s_\beta`. This target distribution is called *local*
Boltzmann since :math:`W` is only modified locally, where sampling has
taken place. We see that when :math:`s_\beta \approx 0` the histogram
essentially does not grow and the size of the free energy update will
stay at a constant value (as in the original formulation of
metadynamics). Thus, the free energy estimate will not converge, but
continue to fluctuate around the correct value. This illustrates the
inherent coupling between the convergence and choice of target
distribution for this special choice of target. Furthermore note that
when using :math:`\rho=\rho_{\mathrm{Boltz,loc}}` there is no initial
stage (section :ref:`awhinitialstage`). The rescaling of the weight
histogram applied in the initial stage is a global operation, which is
incompatible :math:`\rho_{\mathrm{Boltz,loc}}` only depending locally on
the sampling history.

The target distribution can also be modulated by arbitrary
probability weights

.. math:: \rho(\lambda) = \rho_0(\lambda) w_{\mathrm{user}}(\lambda).
          :label: eqnawhpropweigth

where :math:`w_{\mathrm{user}}(\lambda)` is provided by user data and
in principle :math:`\rho_0(\lambda)` can be any of the target
distributions mentioned above.

Lastly, it is possible to automatically scale the target distribution
(:math:`\rho_0(\lambda)`) based on the AWH friction metric (see
section :ref:`awhfriction`). This implies scaling the target
distribution by the square root of the friction metric
(see :eq:`Eq. %s <eqawhsqrtmetric>`),

.. math:: \rho(\lambda) = \rho_0(\lambda) w_{\mathrm{user}}(\lambda) \sqrt{\det\eta_{\mu\nu}(\lambda)},
          :label: eqnawhmetricopt

where :math:`w_{\mathrm{user}}(\lambda)` can be uniform and
\sqrt{\det\eta_{\mu\nu}(\lambda)} is the square root of the friction metric.
This scaling of the target distribution, increasing the
relative sampling of regions with slower diffusion, should generally lower the
statistical error of the estimated free energy landscape.

This modification is only applied after leaving the initial stage
(section :ref:`awhinitialstage`), if applicable, and is performed when updating
the target distribution, typically when also updating the free energy. The scaling
is based on the relative difference of the local friction metric compared to the
average friction metric (of points that have a non-zero friction metric).

If any histograms have not been sampled enough to have a friction metric they
will be scaled by the average friction metric, i.e., practically unscaled.
Insufficient sampling can result in a too low, but still non-zero, friction metric.
To address that, the scaling down of the target distribution (relative scaling < 1)
is based on the local sampling of each point, so that the target distribution of
points that have not been sampled much yet will be almost unaffected.
Furthermore, the scaling can be limited by a maximum scaling factor
(:mdp:`awh1-target-metric-scaling-limit`). The lower limit
of the scaling is the inverse of the maximum scaling factor.

More information about this scaling can be found in :ref:`194 <reflundborg2023>`.

Scaling the target distribution based on the friction metric
can be combined with Boltzmann or Local-Boltzmann target distributions.
However, this is generally not recommended, due to the risk of feedback loops
between the two adaptive update mechanisms.

Multiple independent or sharing biases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple independent bias potentials may be applied within one
simulation. This only makes sense if the biased coordinates
:math:`\xi^{(1)}`, :math:`\xi^{(2)}`, :math:`\ldots` evolve essentially
independently from one another. A typical example of this would be when
applying an independent bias to each monomer of a protein. Furthermore,
multiple AWH simulations can be launched in parallel, each with a (set
of) independent biases.

If the defined sampling interval is large relative to the diffusion time
of the reaction coordinate, traversing the sampling interval multiple
times as is required by the initial stage
(section :ref:`awhinitialstage`) may take an infeasible mount of
simulation time. In these cases it could be advantageous to parallelize
the work and have a group of multiple “walkers” :math:`\xi^{(i)}(t)`
share a single bias potential. This can be achieved by collecting
samples from all :math:`\xi^{(i)}` of the same sharing group into a
single histogram and update a common free energy estimate. Samples can
be shared between walkers within the simulation and/or between multiple
simulations. However, currently only sharing between simulations is
supported in the code while all biases within a simulation are
independent.

Note that when attempting to shorten the simulation time by using
bias-sharing walkers, care must be taken to ensure the simulations are
still long enough to properly explore and equilibrate all regions of the
sampling interval. To begin, the walkers in a group should be
decorrelated and distributed approximately according to the target
distribution before starting to refine the free energy. This can be
achieved e.g. by “equilibrating” the shared weight histogram before
letting it grow; for instance, :math:`W(\lambda)/N\approx \rho(\lambda)`
with some tolerance.

Furthermore, the “covering” or transition criterion of the initial stage
should to be generalized to detect when the sampling interval has been
collectively traversed. One alternative is to just use the same
criterion as for a single walker (but now with more samples), see
:eq:`Eq. %s <eqawhcovering>`. However, in contrast to the
single walker case this does not ensure that any real transitions across
the sampling interval has taken place; in principle all walkers could be
sampling only very locally and still cover the whole interval. Just as
with a standard umbrella sampling procedure, the free energy may appear
to be converged while in reality simulations sampling closeby
:math:`\lambda` values are sampling disconnected regions of phase space.
A stricter criterion, which helps avoid such issues, is to require that
before a simulation marks a point :math:`\lambda_{\mu}` along dimension
:math:`\mu` as visited, and shares this with the other walkers, also all
points within a certain diameter :math:`D_{\mathrm{cover}}` should have
been visited (i.e. fulfill :eq:`Eq. %s <eqawhcovering>`).
Increasing :math:`D_{\mathrm{cover}}` increases robustness, but may slow
down convergence. For the maximum value of :math:`D_{\mathrm{cover}}`,
equal to the length of the sampling interval, the sampling interval is
considered covered when at least one walker has independently traversed
the sampling interval.

In practice biases are shared by setting :mdp:`awh-share-multisim` to true
and :mdp:`awh1-share-group` (for bias 1) to a non-zero value. Here, bias 1
will be shared between simulations that have the same share group value.
Sharing can be different for bias 1, 2, etc. (although there are
few use cases where this is useful). Technically there are no restrictions
on sharing, apart from that biases that are shared need to have the same
number of grid points and the update intervals should match.
Note that biases can not be shared within a simulation.
The latter could be useful, especially for multimeric proteins, but this
is more difficult to implement.

.. _awhreweight:

Reweighting and combining biased data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often one may want to, post-simulation, calculate the unbiased PMF
:math:`\Phi(u)` of another variable :math:`u(x)`. :math:`\Phi(u)` can be
estimated using :math:`\xi`-biased data by reweighting (“unbiasing”) the
trajectory using the bias potential :math:`U_{n(t)}`, see
:eq:`Eq. %s <eqawhbiaspotential>`. Essentially, one bins the
biased data along :math:`u` and removes the effect of :math:`U_{n(t)}`
by dividing the weight of samples :math:`u(t)` by
:math:`e^{-U_{n(t)}(\xi(t))}`,

.. math:: \hat{\Phi}(u)  = -\ln
          \sum_t 1_u(u(t))e^{U_{n(t)}(\xi(t)} \mathcal{Z}_{n(t)}.
          :label: eqawhunbias

Here the indicator function :math:`1_u` denotes the binning procedure:
:math:`1_u(u') = 1` if :math:`u'` falls into the bin labeled by
:math:`u` and :math:`0` otherwise. The normalization factor
:math:`\mathcal{Z}_n = \int e^{-\Phi(\xi) - U_{n}(\xi)}d \xi` is the
partition function of the extended ensemble. As can be seen
:math:`\mathcal{Z}_n` depends on :math:`\Phi(\xi)`, the PMF of the
(biased) reaction coordinate :math:`\xi` (which is calculated and
written to file by the AWH simulation). It is advisable to use only
final stage data in the reweighting procedure due to the rapid change of
the bias potential during the initial stage. If one would include
initial stage data, one should use the sample weights that are inferred
by the repeated rescaling of the histogram in the initial stage, for the
sake of consistency. Initial stage samples would then in any case be
heavily scaled down relative to final stage samples. Note that
:eq:`Eq. %s <eqawhunbias>` can also be used to combine data
from multiple simulations (by adding another sum also over the
trajectory set). Furthermore, when multiple independent AWH biases have
generated a set of PMF estimates :math:`\{\hat{\Phi}^{(i)}(\xi)\}`, a
combined best estimate :math:`\hat{\Phi}(\xi)` can be obtained by
applying self-consistent exponential averaging. More details on this
procedure and a derivation of :eq:`Eq. %s <eqawhunbias>`
(using slightly different notation) can be found in :ref:`143 <reflindahl2017sequence>`.

.. _awhfriction:

The friction metric
~~~~~~~~~~~~~~~~~~~

During the AWH simulation, the following time-integrated force
correlation function is calculated,

.. math:: \eta_{\mu\nu}(\lambda) =
          \beta
          \int_0^\infty
          \frac{
          \left<{\delta \mathcal{F}_{\mu}(x(t),\lambda)
          \delta \mathcal{F}_\nu(x(0),\lambda)
          \omega(\lambda|x(t)) \omega(\lambda|x(0))}\right>}
          {\left<{\omega^2(\lambda | x)}\right>}
          dt.
          :label: eqawhmetric

Here
:math:`\mathcal F_\mu(x,\lambda) = k_\mu (\xi_\mu(x) - \lambda_\mu)` is
the force along dimension :math:`\mu` from an harmonic potential
centered at :math:`\lambda` and
:math:`\delta \mathcal F_{\mu}(x,\lambda) = \mathcal F_{\mu}(x,\lambda) - \left<{\mathcal F_\mu(x,\lambda)}\right>`
is the deviation of the force. The factors :math:`\omega(\lambda|x(t))`,
see :eq:`Eq %s <eqawhomega>`, reweight the samples.
:math:`\eta_{\mu\nu}(\lambda)` is a friction
tensor :ref:`186 <reflindahl2018>` and :ref:`144 <refsivak2012thermodynamic>`.
The diffusion matrix, on the flattened landscape, is equal to :math:`k_B T` times
the inverse of the friction metrix tensor:

.. math:: \mathbf{D}(\lambda) = k_B T \mathbf{\eta}^{-1}(\lambda).
          :label: eqawhdiffusionmatrix

A measure of sampling (in)efficiency at each
:math:`\lambda` is given by

.. math:: \eta^{\frac{1}{2}}(\lambda) = \sqrt{\det\eta_{\mu\nu}(\lambda)}.
          :label: eqawhsqrtmetric

A large value of :math:`\eta^{\frac{1}{2}}(\lambda)` indicates slow
dynamics and long correlation times, which may require more sampling.

.. _awhusage:

Limitations
~~~~~~~~~~~

The only real limitation of the AWH implementation, apart from the not uncommon
practical issue that the method might not converge sufficiently fast, is a limit
on the maximum free energy difference. This limit is set to :math:`700 k_B T`, because
:math:`e^700` is close to the maximum value that can be accurately represented by
a double-precision floating-point value. For physical reaction coordinates, this
is not a limit in practice. This does limit the range of applications for
alchemical coordinates. For instance, hydration free-energies of divalent
cations with a pair of monovalent anions can exceed this limit. The limit can
also be exceeded when decoupling large molecules from solvent, but this often
coincides with the limit where the sampling becomes problematic.

Usage
~~~~~

AWH stores data in the energy file (:ref:`edr`) with a frequency set by the
user. The data – the PMF, the convolved bias, distributions of the
:math:`\lambda` and :math:`\xi` coordinates, etc. – can be extracted
after the simulation using the :ref:`gmx awh` tool. Furthermore, the trajectory
of the reaction coordinate :math:`\xi(t)` is printed to the pull output
file :math:`{\tt pullx.xvg}`. The log file (:ref:`log`) also contains
information; check for messages starting with “awh”, they will tell you
about covering and potential sampling issues.

Setting the initial update size
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The initial value of the weight histogram size :math:`N` sets the
initial update size (and the rate of change of the bias). When :math:`N`
is kept constant, like in the initial stage, the average variance of the
free energy scales as :math:`\varepsilon^2 \sim 1/(ND)`
:ref:`137 <reflindahl2014accelerated>`, for a simple model system with constant diffusion
:math:`D` along the reaction coordinate. This provides a ballpark
estimate used by AWH to initialize :math:`N` in terms of more meaningful
quantities

.. math:: \frac{1}{N_0} = \frac{1}{N_0(\varepsilon_0, D)} = \Delta
	  t_\mathrm{sample} \max_d \frac{2D_d}{L_d^2} \varepsilon_0^2
          :label: eqawhn0

where :math:`L_d` is the length of the interval and :math:`D_d` is
the diffusion constant along dimension :math:`d` of the AWH bias.
For one dimension, :math:`L^2/2D` is the average time to diffuse
over a distance of :math:`L`. We then takes the maximum crossing
time over all dimensions involved in the bias.
Essentially, this formula tells us that a slower system (small :math:`D`)
requires more samples (larger :math:`N^0`) to attain the same level of
accuracy (:math:`\varepsilon_0`) at a given sampling rate. Conversely,
for a system of given diffusion, how to choose the initial biasing rate
depends on how good the initial accuracy is. Both the initial error
:math:`\varepsilon_0` and the diffusion :math:`D` only need to be
roughly estimated or guessed. In the typical case, one would only tweak
the :math:`D` parameter, and use a default value for
:math:`\varepsilon_0`. For good convergence, :math:`D` should be chosen
as large as possible (while maintaining a stable system) giving large
initial bias updates and fast initial transitions. Choosing :math:`D`
too small can lead to slow initial convergence. It may be a good idea to
run a short trial simulation and after the first covering check the
maximum free energy difference of the PMF estimate. If this is much
larger than the expected magnitude of the free energy barriers that
should be crossed, then the system is probably being pulled too hard and
:math:`D` should be decreased. An accurate estimate of the diffusion
can be obtained from an AWH simulation with the :ref:`gmx awh` tool.
:math:`\varepsilon_0` on the other hand, should be a rough estimate
of the initial error.

Estimating errors
^^^^^^^^^^^^^^^^^

As with any adaptive method, estimating errors for AWH is difficult from
data of a single simulation only. We are looking into methods to do this.
For now, the only safe way to estimate errors is to run multiple completely
independent simulations and compute a standard error estimate.
Note that for the simulations to be really independent, they should start
from different, equilibrated states along the reaction coordinate(s).
In practice, this is often difficult to achieve, in particular in the common
case that you only know the starting state along the reaction coordinate.
The exit from the initial phase of AWH is designed such that, in most cases,
such systematic errors are as small as the noise when exiting the initial phase,
but it cannot be excluded that some effects are still present.

Tips for efficient sampling
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The force constant :math:`k` should be larger than the curvature of the
PMF landscape. If this is not the case, the distributions of the
reaction coordinate :math:`\xi` and the reference coordinate
:math:`\lambda`, will differ significantly and warnings will be printed
in the log file. One can choose :math:`k` as large as the time step
supports. This will necessarily increase the number of points of the
discretized sampling interval :math:`I`. In general however, it should
not affect the performance of the simulation noticeably because the AWH
update is implemented such that only sampled points are accessed at free
energy update time.

For an alchemical free-energy dimension, AWH accesses all
:math:`\lambda` points at every sampling step. Because the number of
:math:`\lambda` points is usually far below 100, there is no significant
cost to this in the AWH method itself. However, foreign energy differences
need to be computed for every :math:`\lambda` value used, which can
become somewhat costly.

As with any method, the choice of reaction coordinate(s) is critical. If
a single reaction coordinate does not suffice, identifying a second
reaction coordinate and sampling the two-dimensional landscape may help.
In this case, using a target distribution with a free energy cutoff (see
:eq:`Eq. %s <eqawhrhocut>`) might be required to avoid
sampling uninteresting regions of very high free energy. Obtaining
accurate free energies for reaction coordinates of much higher
dimensionality than 3 or possibly 4 is generally not feasible.

Monitoring the transition rate of :math:`\xi(t)`, across the sampling
interval is also advisable. For reliable statistics (e.g. when
reweighting the trajectory as described in section :ref:`awhreweight`),
one would generally want to observe at least a few transitions after
having exited the initial stage. Furthermore, if the dynamics of the
reaction coordinate suddenly changes, this may be a sign of e.g. a
reaction coordinate problem.

Difficult regions of sampling may also be detected by calculating the
friction tensor :math:`\eta_{\mu\nu}(\lambda)` in the sampling interval,
see section :ref:`awhfriction`. :math:`\eta_{\mu\nu}(\lambda)` as well
as the sampling efficiency measure :math:`\eta^{\frac{1}{2}}(\lambda)`
(:eq:`Eq. %s <eqawhsqrtmetric>`) are written to the energy file and can be
extracted with :ref:`gmx awh`. A high peak in
:math:`\eta^{\frac{1}{2}}(\lambda)` indicates that this region requires
longer time to sample properly.
