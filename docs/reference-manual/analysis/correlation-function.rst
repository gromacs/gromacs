Correlation functions
---------------------

Theory of correlation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The theory of correlation functions is well established \ :ref:`108 <refAllen87>`.
We describe here the implementation of the various
correlation function flavors in the |Gromacs| code. The definition of the
autocorrelation function (ACF) :math:`C_f(t)` for a property
:math:`f(t)` is:

.. math:: C_f(t)  ~=~     \left\langle f(\xi) f(\xi+t)\right\rangle_{\xi}
          :label: eqncorr

where the notation on the right hand side indicates averaging over
:math:`\xi`, *i.e.* over time origins. It is also possible to compute
cross-correlation function from two properties :math:`f(t)` and
:math:`g(t)`:

.. math:: C_{fg}(t) ~=~   \left\langle f(\xi) g(\xi+t)\right\rangle_{\xi}
          :label: eqncrosscorr

however, in |Gromacs| there is no standard mechanism to do this
(**note:** you can use the ``xmgr`` program to compute cross correlations).
The integral of the correlation function over time is the correlation
time :math:`\tau_f`:

.. math:: \tau_f  ~=~     \int_0^{\infty} C_f(t) {\rm d} t
          :label: eqncorrtime

In practice, correlation functions are calculated based on data points
with discrete time intervals :math:`\Delta`\ t, so that the ACF from an
MD simulation is:

.. math:: C_f(j\Delta t)  ~=~     \frac{1}{N-j}\sum_{i=0}^{N-1-j} f(i\Delta t) f((i+j)\Delta t)
          :label: eqncorrmd

where :math:`N` is the number of available time frames for the
calculation. The resulting ACF is obviously only available at time
points with the same interval :math:`\Delta`\ t. Since, for many
applications, it is necessary to know the short time behavior of the ACF
(*e.g.* the first 10 ps) this often means that we have to save the data
with intervals much shorter than the time scale of interest. Another
implication of :eq:`eqn. %s <eqncorrmd>` is that in principle we can not compute
all points of the ACF with the same accuracy, since we have :math:`N-1`
data points for :math:`C_f(\Delta t)` but only 1 for
:math:`C_f((N-1)\Delta t)`. However, if we decide to compute only an ACF
of length :math:`M\Delta t`, where :math:`M \leq N/2` we can compute all
points with the same statistical accuracy:

.. math:: C_f(j\Delta t)  ~=~ \frac{1}{M}\sum_{i=0}^{N-1-M} f(i\Delta t)f((i+j)\Delta t)
          :label: eqncorrstataccuracy

Here of course :math:`j < M`. :math:`M` is sometimes referred to as the
time lag of the correlation function. When we decide to do this, we
intentionally do not use all the available points for very short time
intervals (:math:`j << M`), but it makes it easier to interpret the
results. Another aspect that may not be neglected when computing ACFs
from simulation is that usually the time origins :math:`\xi`
(:eq:`eqn. %s <eqncorr>`) are not statistically independent, which may introduce
a bias in the results. This can be tested using a block-averaging
procedure, where only time origins with a spacing at least the length of
the time lag are included, *e.g.* using :math:`k` time origins with
spacing of :math:`M\Delta t` (where :math:`kM \leq N`):

.. math:: C_f(j\Delta t)  ~=~ \frac{1}{k}\sum_{i=0}^{k-1} f(iM\Delta t)f((iM+j)\Delta t)
          :label: eqncorrblockaveraging

However, one needs very long simulations to get good accuracy this way,
because there are many fewer points that contribute to the ACF.

Using FFT for computation of the ACF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The computational cost for calculating an ACF according to
:eq:`eqn. %s <eqncorrmd>` is proportional to :math:`N^2`, which is considerable.
However, this can be improved by using fast Fourier transforms to do the
convolution \ :ref:`108 <refAllen87>`.

Special forms of the ACF
~~~~~~~~~~~~~~~~~~~~~~~~

There are some important varieties on the ACF, *e.g.* the ACF of a
vector :math:`\mathbf{p}`:

.. math:: C_{\mathbf{p}}(t) ~=~       \int_0^{\infty} P_n(\cos\angle\left(\mathbf{p}(\xi),\mathbf{p}(\xi+t)\right) {\rm d} \xi
          :label: eqncorrleg

where :math:`P_n(x)` is the :math:`n^{th}` order Legendre
polynomial. [1]_ Such correlation times can actually be obtained
experimentally using *e.g.* NMR or other relaxation experiments. |Gromacs|
can compute correlations using the 1\ :math:`^{st}` and 2\ :math:`^{nd}`
order Legendre polynomial (:eq:`eqn. %s <eqncorrleg>`). This can also be used
for rotational autocorrelation (:ref:`gmx rotacf`) and dipole autocorrelation
(:ref:`gmx dipoles <gmx dipoles>`).

In order to study torsion angle dynamics, we define a dihedral
autocorrelation function as \ :ref:`159 <refSpoel97a>`:

.. math:: C(t)    ~=~     \left\langle \cos(\theta(\tau)-\theta(\tau+t))\right\rangle_{\tau}
          :label: eqncoenk

**Note** that this is not a product of two functions as is generally
used for correlation functions, but it may be rewritten as the sum of
two products:

.. math:: C(t)    ~=~     \left\langle\cos(\theta(\tau))\cos(\theta(\tau+t))\,+\,\sin(\theta(\tau))\sin(\theta(\tau+t))\right\rangle_{\tau}
          :label: eqncot

Some Applications
~~~~~~~~~~~~~~~~~

The program :ref:`gmx velacc <gmx velacc>`
calculates the *velocity autocorrelation function*.

.. math:: C_{\mathbf{v}} (\tau) ~=~ \langle {\mathbf{v}}_i(\tau) \cdot {\mathbf{v}}_i(0) \rangle_{i \in A}
          :label: eqnvelocityautocorr

The self diffusion coefficient can be calculated using the Green-Kubo
relation \ :ref:`108 <refAllen87>`:

.. math:: D_A ~=~ {1\over 3} \int_0^{\infty} \langle {\bf v}_i(t) \cdot {\bf v}_i(0) \rangle_{i \in A} \; dt
          :label: eqndiffcoeff

which is just the integral of the velocity autocorrelation function.
There is a widely-held belief that the velocity ACF converges faster
than the mean square displacement (sec. :ref:`msd`), which can also be
used for the computation of diffusion constants. However, Allen &
Tildesley \ :ref:`108 <refAllen87>` warn us that the long-time
contribution to the velocity ACF can not be ignored, so care must be
taken.

Another important quantity is the dipole correlation time. The *dipole
correlation function* for particles of type :math:`A` is calculated as
follows by :ref:`gmx dipoles <gmx dipoles>`:

.. math:: C_{\mu} (\tau) ~=~
          \langle {\bf \mu}_i(\tau) \cdot {\bf \mu}_i(0) \rangle_{i \in A}
          :label: eqndipolecorrfunc

with :math:`{\bf \mu}_i = \sum_{j \in i} {\bf r}_j q_j`. The dipole
correlation time can be computed using :eq:`eqn. %s <eqncorrtime>`. 
For some applications
see (**???**).

The viscosity of a liquid can be related to the correlation time of the
Pressure tensor
:math:`\mathbf{P}` :ref:`160 <refPSmith93c>`,
:ref:`161 <refBalasubramanian96>`. :ref:`gmx energy` can compute the viscosity,
but this is not very accurate \ :ref:`149 <refHess2002a>`, and actually
the values do not converge.


.. [1]
   :math:`P_0(x) = 1`, :math:`P_1(x) = x`, :math:`P_2(x) = (3x^2-1)/2`

