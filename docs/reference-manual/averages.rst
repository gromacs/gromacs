Averages and fluctuations
=========================

Formulae for averaging
----------------------

**Note:** this section was taken from ref \ :ref:`179 <refGunsteren94a>`.

When analyzing a MD trajectory averages :math:`\left<x\right>` and
fluctuations

.. math::  \left<(\Delta x)^2\right>^{{\frac{1}{2}}} ~=~ \left<[x-\left<x\right>]^2\right>^{{\frac{1}{2}}}
           :label: eqnvar0

of a quantity :math:`x` are to be computed. The variance
:math:`\sigma_x` of a series of N\ :math:`_x` values, {:math:`x_i`}, can
be computed from

.. math:: \sigma_x~=~ \sum_{i=1}^{N_x} x_i^2 ~-~  \frac{1}{N_x}\left(\sum_{i=1}^{N_x}x_i\right)^2
          :label: eqnvar1

Unfortunately this formula is numerically not very accurate, especially
when :math:`\sigma_x^{{\frac{1}{2}}}` is small compared to the values of
:math:`x_i`. The following (equivalent) expression is numerically more
accurate

.. math:: \sigma_x ~=~ \sum_{i=1}^{N_x} [x_i  - \left<x\right>]^2
          :label: eqnvar1equivalent

with

.. math:: \left<x\right> ~=~ \frac{1}{N_x} \sum_{i=1}^{N_x} x_i
          :label: eqnvar2

Using :eq:`eqns. %s <eqnvar1>` and
:eq:`%s <eqnvar2>` one has to go through the series of
:math:`x_i` values twice, once to determine :math:`\left<x\right>` and
again to compute :math:`\sigma_x`, whereas
:eq:`eqn. %s <eqnvar0>` requires only one sequential scan of
the series {:math:`x_i`}. However, one may cast
:eq:`eqn. %s <eqnvar1>` in another form, containing partial
sums, which allows for a sequential update algorithm. Define the partial
sum

.. math:: X_{n,m} ~=~ \sum_{i=n}^{m} x_i
          :label: eqnpartialsum

and the partial variance

.. math::  \sigma_{n,m} ~=~ \sum_{i=n}^{m}  \left[x_i - \frac{X_{n,m}}{m-n+1}\right]^2  
           :label: eqnsigma

It can be shown that

.. math::  X_{n,m+k} ~=~  X_{n,m} + X_{m+1,m+k}         
           :label: eqnXpartial

and

.. math:: \begin{aligned}
          \sigma_{n,m+k} &=& \sigma_{n,m} + \sigma_{m+1,m+k} + \left[~\frac {X_{n,m}}{m-n+1} - \frac{X_{n,m+k}}{m+k-n+1}~\right]^2~* \nonumber\\
          && ~\frac{(m-n+1)(m+k-n+1)}{k}
          \end{aligned}
          :label: eqnvarpartial

For :math:`n=1` one finds

.. math:: \sigma_{1,m+k} ~=~ \sigma_{1,m} + \sigma_{m+1,m+k}~+~
          \left[~\frac{X_{1,m}}{m} - \frac{X_{1,m+k}}{m+k}~\right]^2~ \frac{m(m+k)}{k}
          :label: eqnsig1

and for :math:`n=1` and :math:`k=1`
:eq:`eqn. %s <eqnvarpartial>` becomes

.. math:: \begin{aligned}
          \sigma_{1,m+1}  &=& \sigma_{1,m} + 
          \left[\frac{X_{1,m}}{m} - \frac{X_{1,m+1}}{m+1}\right]^2 m(m+1)\\
          &=& \sigma_{1,m} + 
          \frac {[~X_{1,m} - m x_{m+1}~]^2}{m(m+1)}
          \end{aligned}
          :label: eqnsimplevar0

where we have used the relation

.. math:: X_{1,m+1} ~=~  X_{1,m} + x_{m+1}                       
          :label: eqnsimplevar1

Using formulae :eq:`eqn. %s <eqnsimplevar0>` and
:eq:`eqn. %s <eqnsimplevar1>` the average

.. math:: \left<x\right> ~=~ \frac{X_{1,N_x}}{N_x}
          :label: eqnfinalaverage

and the fluctuation

.. math:: \left<(\Delta x)^2\right>^{{\frac{1}{2}}} = \left[\frac {\sigma_{1,N_x}}{N_x}\right]^{{\frac{1}{2}}}
          :label: eqnfinalfluctuation

can be obtained by one sweep through the data.

Implementation
--------------

In |Gromacs| the instantaneous energies :math:`E(m)` are stored in the
:ref:`energy file <edr>`, along with the values of :math:`\sigma_{1,m}` and
:math:`X_{1,m}`. Although the steps are counted from 0, for the energy
and fluctuations steps are counted from 1. This means that the equations
presented here are the ones that are implemented. We give somewhat
lengthy derivations in this section to simplify checking of code and
equations later on.

Part of a Simulation
~~~~~~~~~~~~~~~~~~~~

It is not uncommon to perform a simulation where the first part, *e.g.*
100 ps, is taken as equilibration. However, the averages and
fluctuations as printed in the :ref:`log file <log>` are computed over the whole
simulation. The equilibration time, which is now part of the simulation,
may in such a case invalidate the averages and fluctuations, because
these numbers are now dominated by the initial drift towards
equilibrium.

Using :eq:`eqns. %s <eqnXpartial>` and
:eq:`%s <eqnvarpartial>` the average and standard deviation
over part of the trajectory can be computed as:

.. math:: \begin{aligned}
          X_{m+1,m+k}     &=& X_{1,m+k} - X_{1,m}                 \\
          \sigma_{m+1,m+k} &=& \sigma_{1,m+k}-\sigma_{1,m} - \left[~\frac{X_{1,m}}{m} - \frac{X_{1,m+k}}{m+k}~\right]^{2}~ \frac{m(m+k)}{k}\end{aligned}
          :label: eqnaveragesimpart

or, more generally (with :math:`p \geq 1` and :math:`q \geq p`):

.. math:: \begin{aligned}
          X_{p,q}         &=&     X_{1,q} - X_{1,p-1}     \\
          \sigma_{p,q}    &=&     \sigma_{1,q}-\sigma_{1,p-1} - \left[~\frac{X_{1,p-1}}{p-1} - \frac{X_{1,q}}{q}~\right]^{2}~ \frac{(p-1)q}{q-p+1}\end{aligned}
          :label: eqnaveragesimpartgeneral

**Note** that implementation of this is not entirely trivial, since
energies are not stored every time step of the simulation. We therefore
have to construct :math:`X_{1,p-1}` and :math:`\sigma_{1,p-1}` from the
information at time :math:`p` using :eq:`eqns. %s <eqnsimplevar0>` and
:eq:`%s <eqnsimplevar1>`:

.. math:: \begin{aligned}
          X_{1,p-1}       &=&     X_{1,p} - x_p   \\
          \sigma_{1,p-1}  &=&     \sigma_{1,p} -  \frac {[~X_{1,p-1} - (p-1) x_{p}~]^2}{(p-1)p}\end{aligned}
          :label: eqnfinalaveragesimpartnote

Combining two simulations
~~~~~~~~~~~~~~~~~~~~~~~~~

Another frequently occurring problem is, that the fluctuations of two
simulations must be combined. Consider the following example: we have
two simulations (A) of :math:`n` and (B) of :math:`m` steps, in which
the second simulation is a continuation of the first. However, the
second simulation starts numbering from 1 instead of from :math:`n+1`.
For the partial sum this is no problem, we have to add :math:`X_{1,n}^A`
from run A:

.. math::  X_{1,n+m}^{AB} ~=~ X_{1,n}^A + X_{1,m}^B
           :label: eqnpscomb

When we want to compute the partial variance from the two components we
have to make a correction :math:`\Delta\sigma`:

.. math:: \sigma_{1,n+m}^{AB} ~=~ \sigma_{1,n}^A + \sigma_{1,m}^B +\Delta\sigma
          :label: eqnscombcorr

if we define :math:`x_i^{AB}` as the combined and renumbered set of
data points we can write:

.. math:: \sigma_{1,n+m}^{AB} ~=~ \sum_{i=1}^{n+m}  \left[x_i^{AB} - \frac{X_{1,n+m}^{AB}}{n+m}\right]^2
          :label: eqnpscombpoints

and thus

.. math:: \sum_{i=1}^{n+m}  \left[x_i^{AB} - \frac{X_{1,n+m}^{AB}}{n+m}\right]^2  ~=~
          \sum_{i=1}^{n}  \left[x_i^{A} - \frac{X_{1,n}^{A}}{n}\right]^2  +
          \sum_{i=1}^{m}  \left[x_i^{B} - \frac{X_{1,m}^{B}}{m}\right]^2  +\Delta\sigma
          :label: eqnpscombresult

or

.. math:: \begin{aligned}
          \sum_{i=1}^{n+m}  \left[(x_i^{AB})^2 - 2 x_i^{AB}\frac{X^{AB}_{1,n+m}}{n+m} + \left(\frac{X^{AB}_{1,n+m}}{n+m}\right)^2  \right] &-& \nonumber \\
          \sum_{i=1}^{n}  \left[(x_i^{A})^2 - 2 x_i^{A}\frac{X^A_{1,n}}{n} + \left(\frac{X^A_{1,n}}{n}\right)^2  \right] &-& \nonumber \\
          \sum_{i=1}^{m}  \left[(x_i^{B})^2 - 2 x_i^{B}\frac{X^B_{1,m}}{m} + \left(\frac{X^B_{1,m}}{m}\right)^2  \right] &=& \Delta\sigma\end{aligned}
          :label: eqnpscombresult2

all the :math:`x_i^2` terms drop out, and the terms independent of the
summation counter :math:`i` can be simplified:

.. math:: \begin{aligned}
          \frac{\left(X^{AB}_{1,n+m}\right)^2}{n+m} \,-\, 
          \frac{\left(X^A_{1,n}\right)^2}{n} \,-\, 
          \frac{\left(X^B_{1,m}\right)^2}{m} &-& \nonumber \\
          2\,\frac{X^{AB}_{1,n+m}}{n+m}\sum_{i=1}^{n+m}x_i^{AB} \,+\,
          2\,\frac{X^{A}_{1,n}}{n}\sum_{i=1}^{n}x_i^{A} \,+\,
          2\,\frac{X^{B}_{1,m}}{m}\sum_{i=1}^{m}x_i^{B} &=& \Delta\sigma\end{aligned}
          :label: eqnpscombsimp

we recognize the three partial sums on the second line and use
:eq:`eqn. %s <eqnpscomb>` to obtain:

.. math:: \Delta\sigma ~=~ \frac{\left(mX^A_{1,n} - nX^B_{1,m}\right)^2}{nm(n+m)}
          :label: eqnpscombused

if we check this by inserting :math:`m=1` we get back
:eq:`eqn. %s <eqnsimplevar0>`

Summing energy terms
~~~~~~~~~~~~~~~~~~~~

The :ref:`gmx energy <gmx energy>` program
can also sum energy terms into one, *e.g.* potential + kinetic = total.
For the partial averages this is again easy if we have :math:`S` energy
components :math:`s`:

.. math::  X_{m,n}^S ~=~ \sum_{i=m}^n \sum_{s=1}^S x_i^s ~=~ \sum_{s=1}^S \sum_{i=m}^n x_i^s ~=~ \sum_{s=1}^S X_{m,n}^s
           :label: eqnsumterms

For the fluctuations it is less trivial again, considering for example
that the fluctuation in potential and kinetic energy should cancel.
Nevertheless we can try the same approach as before by writing:

.. math:: \sigma_{m,n}^S ~=~ \sum_{s=1}^S \sigma_{m,n}^s + \Delta\sigma
          :label: eqnsigmatermsfluct

if we fill in :eq:`eqn. %s <eqnsigma>`:

.. math:: \sum_{i=m}^n \left[\left(\sum_{s=1}^S x_i^s\right) - \frac{X_{m,n}^S}{m-n+1}\right]^2 ~=~
          \sum_{s=1}^S \sum_{i=m}^n \left[\left(x_i^s\right) - \frac{X_{m,n}^s}{m-n+1}\right]^2 + \Delta\sigma
          :label: eqnsigmaterms

which we can expand to:

.. math:: \begin{aligned}
          &~&\sum_{i=m}^n \left[\sum_{s=1}^S (x_i^s)^2 + \left(\frac{X_{m,n}^S}{m-n+1}\right)^2 -2\left(\frac{X_{m,n}^S}{m-n+1}\sum_{s=1}^S x_i^s + \sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'} \right)\right]    \nonumber \\
          &-&\sum_{s=1}^S \sum_{i=m}^n \left[(x_i^s)^2 - 2\,\frac{X_{m,n}^s}{m-n+1}\,x_i^s + \left(\frac{X_{m,n}^s}{m-n+1}\right)^2\right] ~=~\Delta\sigma \end{aligned}
          :label: eqnsimtermsexpanded

the terms with :math:`(x_i^s)^2` cancel, so that we can simplify to:

.. math:: \begin{aligned}
          &~&\frac{\left(X_{m,n}^S\right)^2}{m-n+1} -2 \frac{X_{m,n}^S}{m-n+1}\sum_{i=m}^n\sum_{s=1}^S x_i^s -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, -        \nonumber \\
          &~&\sum_{s=1}^S \sum_{i=m}^n \left[- 2\,\frac{X_{m,n}^s}{m-n+1}\,x_i^s + \left(\frac{X_{m,n}^s}{m-n+1}\right)^2\right] ~=~\Delta\sigma \end{aligned}
          :label: eqnsigmatermssimplefied

or

.. math:: -\frac{\left(X_{m,n}^S\right)^2}{m-n+1}  -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, +  \sum_{s=1}^S \frac{\left(X_{m,n}^s\right)^2}{m-n+1}  ~=~\Delta\sigma
           :label: eqnsigmatermsalternative

If we now expand the first term using
:eq:`eqn. %s <eqnsumterms>` we obtain:

.. math:: -\frac{\left(\sum_{s=1}^SX_{m,n}^s\right)^2}{m-n+1}  -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, +      \sum_{s=1}^S \frac{\left(X_{m,n}^s\right)^2}{m-n+1}  ~=~\Delta\sigma
          :label: eqnsigmatermsfirstexpand

which we can reformulate to:

.. math:: -2\left[\sum_{s=1}^S \sum_{s'=s+1}^S X_{m,n}^s X_{m,n}^{s'}\,+\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\right] ~=~\Delta\sigma
          :label: eqnsigmatermsreformed

or

.. math:: -2\left[\sum_{s=1}^S X_{m,n}^s \sum_{s'=s+1}^S X_{m,n}^{s'}\,+\,\sum_{s=1}^S \sum_{i=m}^nx_i^s \sum_{s'=s+1}^S x_i^{s'}\right] ~=~\Delta\sigma
          :label: eqnsigmatermsreformedalternative

which gives

.. math:: -2\sum_{s=1}^S \left[X_{m,n}^s \sum_{s'=s+1}^S \sum_{i=m}^n x_i^{s'}\,+\,\sum_{i=m}^n x_i^s \sum_{s'=s+1}^S x_i^{s'}\right] ~=~\Delta\sigma
          :label: eqnsigmatermsfinal

Since we need all data points :math:`i` to evaluate this, in general
this is not possible. We can then make an estimate of
:math:`\sigma_{m,n}^S` using only the data points that are available
using the left hand side of :eq:`eqn. %s <eqnsigmaterms>`.
While the average can be computed using all time steps in the
simulation, the accuracy of the fluctuations is thus limited by the
frequency with which energies are saved. Since this can be easily done
with a program such as ``xmgr`` this is not
built-in in |Gromacs|.

.. raw:: latex

    \clearpage


