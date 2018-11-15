Expanded Ensemble
-----------------

In an expanded ensemble simulation \ :ref:`68 <refLyubartsev1992>`, both the
coordinates and the thermodynamic ensemble are treated as configuration
variables that can be sampled over. The probability of any given state
can be written as:

.. math:: P(\vec{x},k) \propto \exp\left(-\beta_k U_k + g_k\right),
          :label: eqnexpandensemble

where :math:`\beta_k = \frac{1}{k_B T_k}` is the :math:`\beta`
corresponding to the :math:`k`\ th thermodynamic state, and :math:`g_k`
is a user-specified weight factor corresponding to the :math:`k`\ th
state. This space is therefore a *mixed*, *generalized*, or *expanded*
ensemble which samples from multiple thermodynamic ensembles
simultaneously. :math:`g_k` is chosen to give a specific weighting of
each subensemble in the expanded ensemble, and can either be fixed, or
determined by an iterative procedure. The set of :math:`g_k` is
frequently chosen to give each thermodynamic ensemble equal probability,
in which case :math:`g_k` is equal to the free energy in non-dimensional
units, but they can be set to arbitrary values as desired. Several
different algorithms can be used to equilibrate these weights, described
in the mdp option listings.

In |Gromacs|, this space is sampled by alternating sampling in the
:math:`k` and :math:`\vec{x}` directions. Sampling in the
:math:`\vec{x}` direction is done by standard molecular dynamics
sampling; sampling between the different thermodynamics states is done
by Monte Carlo, with several different Monte Carlo moves supported. The
:math:`k` states can be defined by different temperatures, or choices of
the free energy :math:`\lambda` variable, or both. Expanded ensemble
simulations thus represent a serialization of the replica exchange
formalism, allowing a single simulation to explore many thermodynamic
states.
