Replica exchange
----------------

Replica exchange molecular dynamics (REMD) is a method that can be used
to speed up the sampling of any type of simulation, especially if
conformations are separated by relatively high energy barriers. It
involves simulating multiple replicas of the same system at different
temperatures and randomly exchanging the complete state of two replicas
at regular intervals with the probability:

.. math:: P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
          \left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)(U_1 - U_2)
          \right] \right)
          :label: eqnREX

where :math:`T_1` and :math:`T_2` are the reference temperatures and
:math:`U_1` and :math:`U_2` are the instantaneous potential energies of
replicas 1 and 2 respectively. After exchange the velocities are scaled
by :math:`(T_1/T_2)^{\pm0.5}` and a neighbor search is performed the
next step. This combines the fast sampling and frequent barrier-crossing
of the highest temperature with correct Boltzmann sampling at all the
different temperatures \ :ref:`60 <refHukushima96a>`,
:ref:`61 <refSugita99>`. We only attempt exchanges for neighboring temperatures as the
probability decreases very rapidly with the temperature difference. One
should not attempt exchanges for all possible pairs in one step. If, for
instance, replicas 1 and 2 would exchange, the chance of exchange for
replicas 2 and 3 not only depends on the energies of replicas 2 and 3,
but also on the energy of replica 1. In |Gromacs| this is solved by
attempting exchange for all *odd* pairs on *odd* attempts and for all
*even* pairs on *even* attempts. If we have four replicas: 0, 1, 2 and
3, ordered in temperature and we attempt exchange every 1000 steps,
pairs 0-1 and 2-3 will be tried at steps 1000, 3000 etc. and pair 1-2 at
steps 2000, 4000 etc.

How should one choose the temperatures? The energy difference can be
written as:

.. math:: U_1 - U_2 =  N_{df} \frac{c}{2} k_B (T_1 - T_2)
          :label: eqnREXEdiff

where :math:`N_{df}` is the total number of degrees of freedom of one
replica and :math:`c` is 1 for harmonic potentials and around 2 for
protein/water systems. If :math:`T_2 = (1+\epsilon) T_1` the probability
becomes:

.. math:: P(1 \leftrightarrow 2)
            = \exp\left( -\frac{\epsilon^2 c\,N_{df}}{2 (1+\epsilon)} \right)
          \approx \exp\left(-\epsilon^2 \frac{c}{2} N_{df} \right)
          :label: eqnREXprob

Thus for a probability of :math:`e^{-2}\approx 0.135` one obtains
:math:`\epsilon \approx 2/\sqrt{c\,N_{df}}`. With all bonds constrained
one has :math:`N_{df} \approx 2\, N_{atoms}` and thus for :math:`c` = 2
one should choose :math:`\epsilon` as :math:`1/\sqrt{N_{atoms}}`.
However there is one problem when using pressure coupling. The density
at higher temperatures will decrease, leading to higher energy
\ :ref:`62 <refSeibert2005a>`, which should be taken into account. The |Gromacs| website
features a so-called ``REMD calculator``, that lets you type in the
temperature range and the number of atoms, and based on that proposes a
set of temperatures.

An extension to the REMD for the isobaric-isothermal ensemble was
proposed by Okabe et al. :ref:`63 <refOkabe2001a>`. In this work the
exchange probability is modified to:

.. math:: P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
          \left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)(U_1 - U_2) +
          \left(\frac{P_1}{k_B T_1} - \frac{P_2}{k_B T_2}\right)\left(V_1-V_2\right)
          \right] \right)
          :label: eqnREXexchangeprob

where :math:`P_1` and :math:`P_2` are the respective reference
pressures and :math:`V_1` and :math:`V_2` are the respective
instantaneous volumes in the simulations. In most cases the differences
in volume are so small that the second term is negligible. It only plays
a role when the difference between :math:`P_1` and :math:`P_2` is large
or in phase transitions.

Hamiltonian replica exchange is also supported in |Gromacs|. In
Hamiltonian replica exchange, each replica has a different Hamiltonian,
defined by the free energy pathway specified for the simulation. The
exchange probability to maintain the correct ensemble probabilities is:

.. math:: P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
          \left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)((U_1(x_2) - U_1(x_1)) + (U_2(x_1) - U_2(x_2)))
          \right]\right)
          :label: eqnREXcorrectensemble

The separate Hamiltonians are defined by the free energy functionality
of |Gromacs|, with swaps made between the different values of
:math:`\lambda` defined in the mdp file.

Hamiltonian and temperature replica exchange can also be performed
simultaneously, using the acceptance criteria:

.. math:: P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
          \left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)(\frac{U_1(x_2) - U_1(x_1)}{k_B T_1} + \frac{U_2(x_1) - U_2(x_2)}{k_B T_2})
          \right] \right)
          :label: eqnREXacceptance

Gibbs sampling replica exchange has also been implemented in
|Gromacs| :ref:`64 <refChodera2011>`. In Gibbs sampling replica exchange,
all possible pairs are tested for exchange, allowing swaps between
replicas that are not neighbors.

Gibbs sampling replica exchange requires no additional potential energy
calculations. However there is an additional communication cost in Gibbs
sampling replica exchange, as for some permutations, more than one round
of swaps must take place. In some cases, this extra communication cost
might affect the efficiency.

All replica exchange variants are options of the :ref:`mdrun <gmx mdrun>` program. It will
only work when MPI is installed, due to the inherent parallelism in the
algorithm. For efficiency each replica can run on a separate rank. See
the manual page of :ref:`mdrun <gmx mdrun>` on how to use these multinode features.
