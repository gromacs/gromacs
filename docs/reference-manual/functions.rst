.. _ff:

Interaction function and force fields
=====================================

To accommodate the potential functions used in some popular force fields
(see :ref:`ff`), |Gromacs| offers a choice of functions, both for
non-bonded interaction and for dihedral interactions. They are described
in the appropriate subsections.

The potential functions can be subdivided into three parts

#. *Non-bonded*: Lennard-Jones or Buckingham, and Coulomb or modified
   Coulomb. The non-bonded interactions are computed on the basis of a
   neighbor list (a list of non-bonded atoms within a certain radius),
   in which exclusions are already removed.

#. *Bonded*: covalent bond-stretching, angle-bending, improper
   dihedrals, and proper dihedrals. These are computed on the basis of
   fixed lists.

#. *Restraints*: position restraints, angle restraints, distance
   restraints, orientation restraints and dihedral restraints, all based
   on fixed lists.

#. *Applied Forces*: externally applied forces, see
   chapter :ref:`special`.

Non-bonded interactions
-----------------------

Non-bonded interactions in |Gromacs| are pair-additive:

.. math:: V(\mathbf{r}_1,\ldots \mathbf{r}_N) = \sum_{i<j}V_{ij}(\mathbf{r}_ij);

.. math:: \mathbf{F}_i = -\sum_j \frac{dV_{ij}(r_{ij})}{dr_{ij}} \frac{\mathbf{r}_ij}{r_{ij}}

Since the potential only depends on the scalar distance, interactions
will be centro-symmetric, i.e. the vectorial partial force on particle
:math:`i` from the pairwise interaction :math:`V_{ij}(r_{ij})` has the
opposite direction of the partial force on particle :math:`j`. For
efficiency reasons, interactions are calculated by loops over
interactions and updating both partial forces rather than summing one
complete nonbonded force at a time. The non-bonded interactions contain
a repulsion term, a dispersion term, and a Coulomb term. The repulsion
and dispersion term are combined in either the Lennard-Jones (or 6-12
interaction), or the Buckingham (or exp-6 potential). In addition,
(partially) charged atoms act through the Coulomb term.

.. _lj:

The Lennard-Jones interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Lennard-Jones potential :math:`V_{LJ}` between two atoms equals:

.. math::

   V_{LJ}({r_{ij}}) =  \frac{C_{ij}^{(12)}}{{r_{ij}}^{12}} -
                           \frac{C_{ij}^{(6)}}{{r_{ij}}^6}

See also :numref:`Fig. %s <fig-lj>` The parameters :math:`C^{(12)}_{ij}` and
:math:`C^{(6)}_{ij}` depend on pairs of *atom types*; consequently they
are taken from a matrix of LJ-parameters. In the Verlet cut-off scheme,
the potential is shifted by a constant such that it is zero at the
cut-off distance.

.. _fig-lj:

.. figure:: plots/f-lj.*
   :width: 8.00000cm

   The Lennard-Jones interaction.

The force derived from this potential is:

.. math:: \mathbf{F}_i(\mathbf{r}_ij) = \left( 12~\frac{C_{ij}^{(12)}}{{r_{ij}}^{13}} -
                                    6~\frac{C_{ij}^{(6)}}{{r_{ij}}^7} \right) {\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}

The LJ potential may also be written in the following form:

.. math:: V_{LJ}(\mathbf{r}_ij) = 4\epsilon_{ij}\left(\left(\frac{\sigma_{ij}} {{r_{ij}}}\right)^{12}
          - \left(\frac{\sigma_{ij}}{{r_{ij}}}\right)^{6} \right)
          :label: eqnsigeps

In constructing the parameter matrix for the non-bonded LJ-parameters,
two types of combination rules can be used within |Gromacs|, only
geometric averages (type 1 in the input section of the force-field
file):

.. math:: \begin{array}{rcl}
          C_{ij}^{(6)}    &=& \left({C_{ii}^{(6)} \, C_{jj}^{(6)}}\right)^{1/2}    \\
          C_{ij}^{(12)}   &=& \left({C_{ii}^{(12)} \, C_{jj}^{(12)}}\right)^{1/2}
          \end{array}
          :label: eqncomb

or, alternatively the Lorentz-Berthelot rules can be used. An
arithmetic average is used to calculate :math:`\sigma_{ij}`, while a
geometric average is used to calculate :math:`\epsilon_{ij}` (type 2):

.. math:: \begin{array}{rcl}
          \sigma_{ij}   &=& \frac{1}{ 2}(\sigma_{ii} + \sigma_{jj})        \\
          \epsilon_{ij} &=& \left({\epsilon_{ii} \, \epsilon_{jj}}\right)^{1/2}
          \end{array}
          :label: eqnlorentzberthelot

finally an geometric average for both parameters can be used (type 3):

.. math:: \begin{array}{rcl}
          \sigma_{ij}   &=& \left({\sigma_{ii} \, \sigma_{jj}}\right)^{1/2}        \\
          \epsilon_{ij} &=& \left({\epsilon_{ii} \, \epsilon_{jj}}\right)^{1/2}
          \end{array}

This last rule is used by the OPLS force field.

Buckingham potential
~~~~~~~~~~~~~~~~~~~~

The Buckingham potential has a more flexible and realistic repulsion
term than the Lennard-Jones interaction, but is also more expensive to
compute. The potential form is:

.. math::

   V_{bh}({r_{ij}}) = A_{ij} \exp(-B_{ij} {r_{ij}}) -
                           \frac{C_{ij}}{{r_{ij}}^6}

.. _fig-bham:

.. figure:: plots/f-bham.*
   :width: 8.00000cm

   The Buckingham interaction.

See also :numref:`Fig. %s <fig-bham>`. The force derived from this is:

.. math::

   \mathbf{F}_i({r_{ij}}) = \left[ A_{ij}B_{ij}\exp(-B_{ij} {r_{ij}}) -
                                    6\frac{C_{ij}}{{r_{ij}}^7} \right] {\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}

.. _coul:

Coulomb interaction
~~~~~~~~~~~~~~~~~~~

The Coulomb interaction between two charge particles is given by:

.. math:: V_c({r_{ij}}) = f \frac{q_i q_j}{{\varepsilon_r}{r_{ij}}}
          :label: eqnvcoul

See also :numref:`Fig. %s <fig-coul>`, where
:math:`f = \frac{1}{4\pi \varepsilon_0} = {138.935\,458}` (see chapter :ref:`defunits`)

.. _fig-coul:

.. figure:: plots/vcrf.*
   :width: 8.00000cm

   The Coulomb interaction (for particles with equal signed charge) with
   and without reaction field. In the latter case
   :math:`{\varepsilon_r}` was 1, :math:`{\varepsilon_{rf}}` was 78, and
   :math:`r_c` was 0.9 nm. The dot-dashed line is the same as the dashed
   line, except for a constant.

The force derived from this potential is:

.. math:: \mathbf{F}_i(\mathbf{r}_ij) = f \frac{q_i q_j}{{\varepsilon_r}{r_{ij}}^2}{\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}

A plain Coulomb interaction should only be used without cut-off or when
all pairs fall within the cut-off, since there is an abrupt, large
change in the force at the cut-off. In case you do want to use a
cut-off, the potential can be shifted by a constant to make the
potential the integral of the force. With the group cut-off scheme, this
shift is only applied to non-excluded pairs. With the Verlet cut-off
scheme, the shift is also applied to excluded pairs and self
interactions, which makes the potential equivalent to a reaction field
with :math:`{\varepsilon_{rf}}=1` (see below).

In |Gromacs| the relative dielectric constant :math:`{\varepsilon_r}` may
be set in the in the input for :ref:`grompp <gmx grompp>`.

.. _coulrf:

Coulomb interaction with reaction field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Coulomb interaction can be modified for homogeneous systems by
assuming a constant dielectric environment beyond the cut-off
:math:`r_c` with a dielectric constant of :math:`{\varepsilon_{rf}}`.
The interaction then reads:

.. math:: V_{crf} ~=~
          f \frac{q_i q_j}{{\varepsilon_r}{r_{ij}}}\left[1+\frac{{\varepsilon_{rf}}-{\varepsilon_r}}{2{\varepsilon_{rf}}+{\varepsilon_r}}
          \,\frac{{r_{ij}}^3}{r_c^3}\right]
          - f\frac{q_i q_j}{{\varepsilon_r}r_c}\,\frac{3{\varepsilon_{rf}}}{2{\varepsilon_{rf}}+{\varepsilon_r}}
          :label: eqnvcrf

in which the constant expression on the right makes the potential zero
at the cut-off :math:`r_c`. For charged cut-off spheres this corresponds
to neutralization with a homogeneous background charge. We can rewrite
:eq:`eqn. %s <eqnvcrf>` for simplicity as

.. math:: V_{crf} ~=~     f \frac{q_i q_j}{{\varepsilon_r}}\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]

with

.. math:: \begin{aligned}
          k_{rf}  &=&     \frac{1}{r_c^3}\,\frac{{\varepsilon_{rf}}-{\varepsilon_r}}{(2{\varepsilon_{rf}}+{\varepsilon_r})}
          \end{aligned}
          :label: eqnkrf

.. math:: \begin{aligned}
          c_{rf}  &=&     \frac{1}{r_c}+k_{rf}\,r_c^2 ~=~ \frac{1}{r_c}\,\frac{3{\varepsilon_{rf}}}{(2{\varepsilon_{rf}}+{\varepsilon_r})}
          \end{aligned}
          :label: eqncrf

For large :math:`{\varepsilon_{rf}}` the :math:`k_{rf}` goes to
:math:`r_c^{-3}/2`, while for :math:`{\varepsilon_{rf}}` =
:math:`{\varepsilon_r}` the correction vanishes. In :numref:`Fig. %s <fig-coul>` the
modified interaction is plotted, and it is clear that the derivative
with respect to :math:`{r_{ij}}` (= -force) goes to zero at the cut-off
distance. The force derived from this potential reads:

.. math:: \mathbf{F}_i(\mathbf{r}_ij) = f \frac{q_i q_j}{{\varepsilon_r}}\left[\frac{1}{{r_{ij}}^2} - 2 k_{rf}{r_{ij}}\right]{\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}
          :label: eqnfcrf

The reaction-field correction should also be applied to all excluded
atoms pairs, including self pairs, in which case the normal Coulomb term
in :eq:`eqns. %s <eqnvcrf>` and :eq:`%s <eqnfcrf>` is absent.

Tironi *et al.* have introduced a generalized reaction field in which
the dielectric continuum beyond the cut-off :math:`r_c` also has an
ionic strength :math:`I` :ref:`71 <refTironi95>`. In this case we can
rewrite the constants :math:`k_{rf}` and :math:`c_{rf}` using the
inverse Debye screening length :math:`\kappa`:

.. math:: \begin{aligned}
          \kappa^2  &=&     
          \frac{2 I \,F^2}{\varepsilon_0 {\varepsilon_{rf}}RT}
          = \frac{F^2}{\varepsilon_0 {\varepsilon_{rf}}RT}\sum_{i=1}^{K} c_i z_i^2     \\
          k_{rf}  &=&     \frac{1}{r_c^3}\,
          \frac{({\varepsilon_{rf}}-{\varepsilon_r})(1 + \kappa r_c) + {\frac{1}{2}}{\varepsilon_{rf}}(\kappa r_c)^2}
          {(2{\varepsilon_{rf}}+ {\varepsilon_r})(1 + \kappa r_c) + {\varepsilon_{rf}}(\kappa r_c)^2}
          \end{aligned}
          :label: eqnkgrf

.. math:: \begin{aligned}
          c_{rf}  &=&     \frac{1}{r_c}\,
          \frac{3{\varepsilon_{rf}}(1 + \kappa r_c + {\frac{1}{2}}(\kappa r_c)^2)}
          {(2{\varepsilon_{rf}}+{\varepsilon_r})(1 + \kappa r_c) + {\varepsilon_{rf}}(\kappa r_c)^2}
          \end{aligned}
          :label: eqncgrf

where :math:`F` is Faraday’s constant, :math:`R` is the ideal gas
constant, :math:`T` the absolute temperature, :math:`c_i` the molar
concentration for species :math:`i` and :math:`z_i` the charge number of
species :math:`i` where we have :math:`K` different species. In the
limit of zero ionic strength (:math:`\kappa=0`) :eq:`eqns. %s <eqnkgrf>` and
:eq:`%s <eqncgrf>` reduce to the simple forms of :eq:`eqns. %s <eqnkrf>` and :eq:`%s <eqncrf>`
respectively.

.. _modnbint:

Modified non-bonded interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In |Gromacs|, the non-bonded potentials can be modified by a shift
function, also called a force-switch function, since it switches the
force to zero at the cut-off. The purpose of this is to replace the
truncated forces by forces that are continuous and have continuous
derivatives at the cut-off radius. With such forces the time integration
produces smaller errors. But note that for Lennard-Jones interactions
these errors are usually smaller than other errors, such as integration
errors at the repulsive part of the potential. For Coulomb interactions
we advise against using a shifted potential and for use of a reaction
field or a proper long-range method such as PME.

There is *no* fundamental difference between a switch function (which
multiplies the potential with a function) and a shift function (which
adds a function to the force or potential) \ :ref:`72 <refSpoel2006a>`. The
switch function is a special case of the shift function, which we apply
to the *force function* :math:`F(r)`, related to the electrostatic or
van der Waals force acting on particle :math:`i` by particle :math:`j`
as:

.. math:: \mathbf{F}_i = c \, F(r_{ij}) \frac{\mathbf{r}_ij}{r_{ij}}

For pure Coulomb or Lennard-Jones interactions
:math:`F(r) = F_\alpha(r) = \alpha \, r^{-(\alpha+1)}`. The switched
force :math:`F_s(r)` can generally be written as:

.. math::

   \begin{array}{rcl}
   F_s(r)~=&~F_\alpha(r)   & r < r_1               \\
   F_s(r)~=&~F_\alpha(r)+S(r)      & r_1 \le r < r_c       \\
   F_s(r)~=&~0             & r_c \le r     
   \end{array}

When :math:`r_1=0` this is a traditional shift function, otherwise it
acts as a switch function. The corresponding shifted potential function
then reads:

.. math:: V_s(r) =  \int^{\infty}_r~F_s(x)\, dx

The |Gromacs| **force switch** function :math:`S_F(r)` should be smooth at
the boundaries, therefore the following boundary conditions are imposed
on the switch function:

.. math::

   \begin{array}{rcl}
   S_F(r_1)          &=&0            \\
   S_F'(r_1)         &=&0            \\
   S_F(r_c)          &=&-F_\alpha(r_c)       \\
   S_F'(r_c)         &=&-F_\alpha'(r_c)
   \end{array}

A 3\ :math:`^{rd}` degree polynomial of the form

.. math:: S_F(r) = A(r-r_1)^2 + B(r-r_1)^3

fulfills these requirements. The constants A and B are given by the
boundary condition at :math:`r_c`:

.. math::

   \begin{array}{rcl}
   A &~=~& -\alpha \, \displaystyle
           \frac{(\alpha+4)r_c~-~(\alpha+1)r_1} {r_c^{\alpha+2}~(r_c-r_1)^2} \\
   B &~=~& \alpha \, \displaystyle
           \frac{(\alpha+3)r_c~-~(\alpha+1)r_1}{r_c^{\alpha+2}~(r_c-r_1)^3}
   \end{array}

Thus the total force function is:

.. math:: F_s(r) = \frac{\alpha}{r^{\alpha+1}} + A(r-r_1)^2 + B(r-r_1)^3

and the potential function reads:

.. math:: V_s(r) = \frac{1}{r^\alpha} - \frac{A}{3} (r-r_1)^3 - \frac{B}{4} (r-r_1)^4 - C

where

.. math:: C =  \frac{1}{r_c^\alpha} - \frac{A}{3} (r_c-r_1)^3 - \frac{B}{4} (r_c-r_1)^4

The |Gromacs| **potential-switch** function :math:`S_V(r)` scales the
potential between :math:`r_1` and :math:`r_c`, and has similar boundary
conditions, intended to produce smoothly-varying potential and forces:

.. math::

   \begin{array}{rcl}
   S_V(r_1)          &=&1 \\
   S_V'(r_1)         &=&0 \\
   S_V''(r_1)        &=&0 \\
   S_V(r_c)          &=&0 \\
   S_V'(r_c)         &=&0 \\
   S_V''(r_c)        &=&0
   \end{array}

The fifth-degree polynomial that has these properties is

.. math:: S_V(r; r_1, r_c) = \frac{1 - 10(r-r_1)^3(r_c-r_1)^2 + 15(r-r_1)^4(r_c-r_1) - 6(r-r_1)}{(r_c-r_1)^5}

This implementation is found in several other simulation
packages,\ :ref:`73 <refOhmine1988>`\ :ref:`75 <refGuenot1993>` but
differs from that in CHARMM.\ :ref:`76 <refSteinbach1994>` Switching the
potential leads to artificially large forces in the switching region,
therefore it is not recommended to switch Coulomb interactions using
this function,\ :ref:`72 <refSpoel2006a>` but switching Lennard-Jones
interactions using this function produces acceptable results.

Modified short-range interactions with Ewald summation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When Ewald summation or particle-mesh Ewald is used to calculate the
long-range interactions, the short-range Coulomb potential must also be
modified. Here the potential is switched to (nearly) zero at the
cut-off, instead of the force. In this case the short range potential is
given by:

.. math:: V(r) = f \frac{\mbox{erfc}(\beta r_{ij})}{r_{ij}} q_i q_j,

where :math:`\beta` is a parameter that determines the relative weight
between the direct space sum and the reciprocal space sum and
erfc\ :math:`(x)` is the complementary error function. For further
details on long-range electrostatics, see sec. :ref:`lrelstat`.

Bonded interactions
-------------------

Bonded interactions are based on a fixed list of atoms. They are not
exclusively pair interactions, but include 3- and 4-body interactions as
well. There are *bond stretching* (2-body), *bond angle* (3-body), and
*dihedral angle* (4-body) interactions. A special type of dihedral
interaction (called *improper dihedral*) is used to force atoms to
remain in a plane or to prevent transition to a configuration of
opposite chirality (a mirror image).

.. _bondpot:

Bond stretching
~~~~~~~~~~~~~~~

.. _harmonicbond:

Harmonic potential
^^^^^^^^^^^^^^^^^^

The bond stretching between two covalently bonded atoms :math:`i` and
:math:`j` is represented by a harmonic potential:

.. _fig-bstretch1:

.. figure:: plots/bstretch.*
   :width: 7.00000cm

   Principle of bond stretching (left), and the bond stretching
   potential (right).

.. math:: V_b~({r_{ij}}) = {\frac{1}{2}}k^b_{ij}({r_{ij}}-b_{ij})^2

See also :numref:`Fig. %s <fig-bstretch1>`, with the force given by:

.. math:: \mathbf{F}_i(\mathbf{r}_ij) = k^b_{ij}({r_{ij}}-b_{ij}) {\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}

.. _g96bond:

Fourth power potential
^^^^^^^^^^^^^^^^^^^^^^

In the GROMOS-96 force field \ :ref:`77 <refgromos96>`, the covalent bond
potential is, for reasons of computational efficiency, written as:

.. math:: V_b~({r_{ij}}) = \frac{1}{4}k^b_{ij}\left({r_{ij}}^2-b_{ij}^2\right)^2

The corresponding force is:

.. math:: \mathbf{F}_i(\mathbf{r}_ij) = k^b_{ij}({r_{ij}}^2-b_{ij}^2)~\mathbf{r}_ij

The force constants for this form of the potential are related to the
usual harmonic force constant :math:`k^{b,\mathrm{harm}}`
(sec. :ref:`bondpot`) as

.. math:: 2 k^b b_{ij}^2 = k^{b,\mathrm{harm}}

The force constants are mostly derived from the harmonic ones used in
GROMOS-87 :ref:`78 <refbiomos>`. Although this form is
computationally more efficient (because no square root has to be
evaluated), it is conceptually more complex. One particular disadvantage
is that since the form is not harmonic, the average energy of a single
bond is not equal to :math:`{\frac{1}{2}}kT` as it is for the normal
harmonic potential.

.. _morsebond:

Morse potential bond stretching
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For some systems that require an anharmonic bond stretching potential,
the Morse potential \ :ref:`79 <refMorse29>` between two atoms *i* and *j* is
available in |Gromacs|. This potential differs from the harmonic potential
in that it has an asymmetric potential well and a zero force at infinite
distance. The functional form is:

.. math:: \displaystyle V_{morse} (r_{ij}) = D_{ij} [1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))]^2,

See also :numref:`Fig. %s <fig-morse>`, and the corresponding force is:

.. math::

   \begin{array}{rcl}
   \displaystyle {\bf F}_{morse} ({\bf r}_{ij})&=&2 D_{ij} \beta_{ij} \exp(-\beta_{ij}(r_{ij}-b_{ij})) * \\
   \displaystyle \: & \: &[1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))] \frac{\displaystyle {\bf r}_{ij}}{\displaystyle r_{ij}},
   \end{array}

where :math:`\displaystyle D_{ij}`  is the depth of the well in
kJ/mol, :math:`\displaystyle \beta_{ij}` defines the steepness of the
well (in nm\ :math:`^{-1}`), and :math:`\displaystyle b_{ij}` is the
equilibrium distance in nm. The steepness parameter
:math:`\displaystyle \beta_{ij}` can be expressed in terms of the reduced mass of the atoms *i* and
*j*, the fundamental vibration frequency :math:`\displaystyle\omega_{ij}` and the well depth :math:`\displaystyle D_{ij}`:

.. math:: \displaystyle \beta_{ij}= \omega_{ij} \sqrt{\frac{\mu_{ij}}{2 D_{ij}}}

and because :math:`\displaystyle \omega = \sqrt{k/\mu}`, one can
rewrite :math:`\displaystyle \beta_{ij}` in terms of the harmonic
force constant :math:`\displaystyle k_{ij}`:

.. math:: \displaystyle \beta_{ij}= \sqrt{\frac{k_{ij}}{2 D_{ij}}}
          :label: eqnbetaij

For small deviations :math:`\displaystyle (r_{ij}-b_{ij})`, one can
approximate the :math:`\displaystyle \exp`-term to first-order using a
Taylor expansion:

.. math:: \displaystyle \exp(-x) \approx 1-x
          :label: eqnexpminx

and substituting :eq:`eqn. %s <eqnbetaij>` and :eq:`eqn. %s <eqnexpminx>` in the
functional form:

.. math::

   \begin{array}{rcl}
   \displaystyle V_{morse} (r_{ij})&=&D_{ij} [1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))]^2\\
   \displaystyle \:&=&D_{ij} [1 - (1 -\sqrt{\frac{k_{ij}}{2 D_{ij}}}(r_{ij}-b_{ij}))]^2\\
   \displaystyle \:&=&\frac{1}{2} k_{ij} (r_{ij}-b_{ij}))^2
   \end{array}

we recover the harmonic bond stretching potential.

.. _fig-morse:

.. figure:: plots/f-morse.*
   :width: 7.00000cm

   The Morse potential well, with bond length 0.15 nm.

Cubic bond stretching potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another anharmonic bond stretching potential that is slightly simpler
than the Morse potential adds a cubic term in the distance to the simple
harmonic form:

.. math:: V_b~({r_{ij}}) = k^b_{ij}({r_{ij}}-b_{ij})^2 + k^b_{ij}k^{cub}_{ij}({r_{ij}}-b_{ij})^3

A flexible water model (based on the SPC water model \ :ref:`80 <refBerendsen81>`)
including a cubic bond stretching potential for the O-H bond was
developed by Ferguson \ :ref:`81 <refFerguson95>`. This model was found to yield a
reasonable infrared spectrum. The Ferguson water model is available in
the |Gromacs| library (``flexwat-ferguson.itp``). It should be noted that the
potential is asymmetric: overstretching leads to infinitely low
energies. The integration timestep is therefore limited to 1 fs.

The force corresponding to this potential is:

.. math:: \mathbf{F}_i(\mathbf{r}_ij) = 2k^b_{ij}({r_{ij}}-b_{ij})~{\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}+ 3k^b_{ij}k^{cub}_{ij}({r_{ij}}-b_{ij})^2~{\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}

FENE bond stretching potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In coarse-grained polymer simulations the beads are often connected by a
FENE (finitely extensible nonlinear elastic) potential \ :ref:`82 <refWarner72>`:

.. math::

   V_{\mbox{\small FENE}}({r_{ij}}) =
     -{\frac{1}{2}}k^b_{ij} b^2_{ij} \log\left(1 - \frac{{r_{ij}}^2}{b^2_{ij}}\right)

The potential looks complicated, but the expression for the force is
simpler:

.. math::

   F_{\mbox{\small FENE}}(\mathbf{r}_ij) =
     -k^b_{ij} \left(1 - \frac{{r_{ij}}^2}{b^2_{ij}}\right)^{-1} \mathbf{r}_ij

At short distances the potential asymptotically goes to a harmonic
potential with force constant :math:`k^b`, while it diverges at distance
:math:`b`.

.. _harmonicangle:

Harmonic angle potential
~~~~~~~~~~~~~~~~~~~~~~~~

The bond-angle vibration between a triplet of atoms :math:`i` -
:math:`j` - :math:`k` is also represented by a harmonic potential on the
angle :math:`{\theta_{ijk}}`

.. _fig-angle:

.. figure:: plots/angle.*
   :width: 7.00000cm

   Principle of angle vibration (left) and the bond angle potential
   (right).

.. math:: V_a({\theta_{ijk}}) = {\frac{1}{2}}k^{\theta}_{ijk}({\theta_{ijk}}-{\theta_{ijk}}^0)^2

As the bond-angle vibration is represented by a harmonic potential, the
form is the same as the bond stretching
(:numref:`Fig. %s <fig-bstretch1>`).

The force equations are given by the chain rule:

.. math::

   \begin{array}{l}
   \mathbf{F}_i    ~=~ -\displaystyle\frac{d V_a({\theta_{ijk}})}{d \mathbf{r}_i}   \\
   \mathbf{F}_k    ~=~ -\displaystyle\frac{d V_a({\theta_{ijk}})}{d \mathbf{r}_k}   \\
   \mathbf{F}_j    ~=~ -\mathbf{F}_i-\mathbf{F}_k
   \end{array}
   ~ \mbox{ ~ where ~ } ~
    {\theta_{ijk}}= \arccos \frac{(\mathbf{r}_ij \cdot \mathbf{r}_{kj})}{r_{ij}r_{kj}}

The numbering :math:`i,j,k` is in sequence of covalently bonded atoms.
Atom :math:`j` is in the middle; atoms :math:`i` and :math:`k` are at
the ends (see :numref:`Fig. %s <fig-angle>`). **Note** that in the input in topology
files, angles are given in degrees and force constants in
kJ/mol/rad\ :math:`^2`.

.. _g96angle:

Cosine based angle potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the GROMOS-96 force field a simplified function is used to represent
angle vibrations:

.. math:: V_a({\theta_{ijk}}) = {\frac{1}{2}}k^{\theta}_{ijk}\left(\cos({\theta_{ijk}}) - \cos({\theta_{ijk}}^0)\right)^2
          :label: eqG96angle

where

.. math:: \cos({\theta_{ijk}}) = \frac{\mathbf{r}_ij\cdot\mathbf{r}_{kj}}{{r_{ij}}r_{kj}}

The corresponding force can be derived by partial differentiation with
respect to the atomic positions. The force constants in this function
are related to the force constants in the harmonic form
:math:`k^{\theta,\mathrm{harm}}` (:ref:`harmonicangle`) by:

.. math:: k^{\theta} \sin^2({\theta_{ijk}}^0) = k^{\theta,\mathrm{harm}}

In the GROMOS-96 manual there is a much more complicated conversion
formula which is temperature dependent. The formulas are equivalent at 0
K and the differences at 300 K are on the order of 0.1 to 0.2%. **Note**
that in the input in topology files, angles are given in degrees and
force constants in kJ/mol.

.. _reb:

Restricted bending potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The restricted bending (ReB) potential \ :ref:`83 <refMonicaGoga2013>` prevents the
bending angle :math:`\theta` from reaching the :math:`180^{\circ}`
value. In this way, the numerical instabilities due to the calculation
of the torsion angle and potential are eliminated when performing
coarse-grained molecular dynamics simulations.

To systematically hinder the bending angles from reaching the
:math:`180^{\circ}` value, the bending potential :eq:`eqn %s <eqG96angle>` is
divided by a :math:`\sin^2\theta` factor:

.. math:: V_{\rm ReB}(\theta_i) = \frac{1}{2} k_{\theta} \frac{(\cos\theta_i - \cos\theta_0)^2}{\sin^2\theta_i}.
          :label: eqReB

:numref:`Figure %s <fig-ReB>` shows the comparison between the ReB potential,
:eq:`%s <eqReB>`, and the standard one :eq:`%s <eqG96angle>`.

.. _fig-ReB:

.. figure:: plots/fig-02.*
   :width: 10.00000cm

   Bending angle potentials: cosine harmonic (solid black line), angle
   harmonic (dashed black line) and restricted bending (red) with the
   same bending constant :math:`k_{\theta}=85` kJ mol\ :math:`^{-1}` and
   equilibrium angle :math:`\theta_0=130^{\circ}`. The orange line
   represents the sum of a cosine harmonic (:math:`k =50` kJ
   mol\ :math:`^{-1}`) with a restricted bending (:math:`k =25` kJ
   mol\ :math:`^{-1}`) potential, both with
   :math:`\theta_0=130^{\circ}`.

The wall of the ReB potential is very repulsive in the region close to
:math:`180^{\circ}` and, as a result, the bending angles are kept within
a safe interval, far from instabilities. The power :math:`2` of
:math:`\sin\theta_i` in the denominator has been chosen to guarantee
this behavior and allows an elegant differentiation:

.. math:: F_{\rm ReB}(\theta_i) = \frac{2k_{\theta}}{\sin^4\theta_i}(\cos\theta_i - \cos\theta_0) (1 - \cos\theta_i\cos\theta_0) \frac{\partial \cos\theta_i}{\partial \vec r_{k}}.
          :label: eqdiffReB

Due to its construction, the restricted bending potential cannot be
used for equilibrium :math:`\theta_0` values too close to
:math:`0^{\circ}` or :math:`180^{\circ}` (from experience, at least
:math:`10^{\circ}` difference is recommended). It is very important
that, in the starting configuration, all the bending angles have to be
in the safe interval to avoid initial instabilities. This bending
potential can be used in combination with any form of torsion potential.
It will always prevent three consecutive particles from becoming
collinear and, as a result, any torsion potential will remain free of
singularities. It can be also added to a standard bending potential to
affect the angle around :math:`180^{\circ}`, but to keep its original
form around the minimum (see the orange curve in :numref:`Fig. %s <fig-ReB>`).

Urey-Bradley potential
~~~~~~~~~~~~~~~~~~~~~~

The Urey-Bradley bond-angle vibration between a triplet of atoms
:math:`i` - :math:`j` - :math:`k` is represented by a harmonic potential
on the angle :math:`{\theta_{ijk}}` and a harmonic correction term on
the distance between the atoms :math:`i` and :math:`k`. Although this
can be easily written as a simple sum of two terms, it is convenient to
have it as a single entry in the topology file and in the output as a
separate energy term. It is used mainly in the CHARMm force
field \ :ref:`84 <refBBrooks83>`. The energy is given by:

.. math:: V_a({\theta_{ijk}}) = {\frac{1}{2}}k^{\theta}_{ijk}({\theta_{ijk}}-{\theta_{ijk}}^0)^2 + {\frac{1}{2}}k^{UB}_{ijk}(r_{ik}-r_{ik}^0)^2

The force equations can be deduced from sections :ref:`harmonicbond`
and :ref:`harmonicangle`.

Bond-Bond cross term
~~~~~~~~~~~~~~~~~~~~

The bond-bond cross term for three particles :math:`i, j, k` forming
bonds :math:`i-j` and :math:`k-j` is given
by \ :ref:`85 <refLawrence2003b>`:

.. math:: V_{rr'} ~=~ k_{rr'} \left(\left|\mathbf{r}_{i}-\mathbf{r}_j\right|-r_{1e}\right) \left(\left|\mathbf{r}_{k}-\mathbf{r}_j\right|-r_{2e}\right)
          :label: eqncrossbb

where :math:`k_{rr'}` is the force constant, and :math:`r_{1e}` and
:math:`r_{2e}` are the equilibrium bond lengths of the :math:`i-j` and
:math:`k-j` bonds respectively. The force associated with this potential
on particle :math:`i` is:

.. math:: \mathbf{F}_{i} = -k_{rr'}\left(\left|\mathbf{r}_{k}-\mathbf{r}_j\right|-r_{2e}\right)\frac{\mathbf{r}_i-\mathbf{r}_j}{\left|\mathbf{r}_{i}-\mathbf{r}_j\right|}

The force on atom :math:`k` can be obtained by swapping :math:`i` and
:math:`k` in the above equation. Finally, the force on atom :math:`j`
follows from the fact that the sum of internal forces should be zero:
:math:`\mathbf{F}_j = -\mathbf{F}_i-\mathbf{F}_k`.

Bond-Angle cross term
~~~~~~~~~~~~~~~~~~~~~

The bond-angle cross term for three particles :math:`i, j, k` forming
bonds :math:`i-j` and :math:`k-j` is given
by \ :ref:`85 <refLawrence2003b>`:

.. math:: V_{r\theta} ~=~ k_{r\theta} \left(\left|\mathbf{r}_{i}-\mathbf{r}_k\right|-r_{3e} \right) \left(\left|\mathbf{r}_{i}-\mathbf{r}_j\right|-r_{1e} + \left|\mathbf{r}_{k}-\mathbf{r}_j\right|-r_{2e}\right)

where :math:`k_{r\theta}` is the force constant, :math:`r_{3e}` is the
:math:`i-k` distance, and the other constants are the same as in
:eq:`Equation %s <eqncrossbb>`. The force associated with the potential on atom
:math:`i` is:

.. math::

   \mathbf{F}_{i} ~=~ -k_{r\theta}
   \left[
   \left(
   \left| \mathbf{r}_{i} - \mathbf{r}_{k}\right|
   -r_{3e}\right)
   \frac{
         \mathbf{r}_{i}-\mathbf{r}_j}
         { \left| \mathbf{r}_{i}-\mathbf{r}_{j}\right| 
         }
   + \left(
     \left| \mathbf{r}_{i}-\mathbf{r}_{j}\right|
   -r_{1e}
   + \left| \mathbf{r}_{k}-\mathbf{r}_{j}\right|
   -r_{2e}\right)
   \frac{
         \mathbf{r}_{i}-\mathbf{r}_{k}}
         {\left| \mathbf{r}_{i}-\mathbf{r}_{k}\right|
         }
   \right]

Quartic angle potential
~~~~~~~~~~~~~~~~~~~~~~~

For special purposes there is an angle potential that uses a fourth
order polynomial:

.. math:: V_q({\theta_{ijk}}) ~=~ \sum_{n=0}^5 C_n ({\theta_{ijk}}-{\theta_{ijk}}^0)^n

.. _imp:

Improper dihedrals
~~~~~~~~~~~~~~~~~~

Improper dihedrals are meant to keep planar groups (*e.g.* aromatic
rings) planar, or to prevent molecules from flipping over to their
mirror images, see :numref:`Fig. %s <fig-imp>`.

.. _fig-imp:

.. figure:: plots/ring-imp.*
        :width: 4.00000cm

        Principle of improper dihedral angles. Out of plane bending for rings.
        The improper dihedral angle :math:`\xi` is defined as the angle between
        planes (i,j,k) and (j,k,l).

.. figure:: plots/subst-im.*
        :width: 3.00000cm

.. figure:: plots/tetra-im.*
        :width: 3.00000cm

        Principle of improper dihedral angles. Out of tetrahedral angle.
        The improper dihedral angle :math:`\xi` is defined
        as the angle between planes (i,j,k) and (j,k,l).

Improper dihedrals: harmonic type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest improper dihedral potential is a harmonic potential; it is
plotted in :numref:`Fig. %s <fig-imps>`.

.. math:: V_{id}(\xi_{ijkl}) = {\frac{1}{2}}k_{\xi}(\xi_{ijkl}-\xi_0)^2

Since the potential is harmonic it is discontinuous, but since the
discontinuity is chosen at 180\ :math:`^\circ` distance from
:math:`\xi_0` this will never cause problems. **Note** that in the input
in topology files, angles are given in degrees and force constants in
kJ/mol/rad\ :math:`^2`.

.. _fig-imps:

.. figure:: plots/f-imps.*
   :width: 10.00000cm

   Improper dihedral potential.

Improper dihedrals: periodic type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This potential is identical to the periodic proper dihedral (see below).
There is a separate dihedral type for this (type 4) only to be able to
distinguish improper from proper dihedrals in the parameter section and
the output.

Proper dihedrals
~~~~~~~~~~~~~~~~

For the normal dihedral interaction there is a choice of either the
GROMOS periodic function or a function based on expansion in powers of
:math:`\cos \phi` (the so-called Ryckaert-Bellemans potential). This
choice has consequences for the inclusion of special interactions
between the first and the fourth atom of the dihedral quadruple. With
the periodic GROMOS potential a special 1-4 LJ-interaction must be
included; with the Ryckaert-Bellemans potential *for alkanes* the 1-4
interactions must be excluded from the non-bonded list. **Note:**
Ryckaert-Bellemans potentials are also used in *e.g.* the OPLS force
field in combination with 1-4 interactions. You should therefore not
modify topologies generated by :ref:`pdb2gmx <gmx pdb2gmx>` in this case.

Proper dihedrals: periodic type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Proper dihedral angles are defined according to the IUPAC/IUB
convention, where :math:`\phi` is the angle between the :math:`ijk` and
the :math:`jkl` planes, with **zero** corresponding to the *cis*
configuration (:math:`i` and :math:`l` on the same side). There are two
dihedral function types in |Gromacs| topology files. There is the standard
type 1 which behaves like any other bonded interactions. For certain
force fields, type 9 is useful. Type 9 allows multiple potential
functions to be applied automatically to a single dihedral in the
``[ dihedral ]`` section when multiple parameters are
defined for the same atomtypes in the ``[ dihedraltypes ]``
section.

.. _fig-pdihf:

.. figure:: plots/f-dih.*
   :width: 7.00000cm

   Principle of proper dihedral angle (left, in *trans* form) and the
   dihedral angle potential (right).

.. math:: V_d(\phi_{ijkl}) = k_{\phi}(1 + \cos(n \phi - \phi_s))

Proper dihedrals: Ryckaert-Bellemans function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| For alkanes, the following proper dihedral potential is often used
  (see :numref:`Fig. %s <fig-rbdih>`):

  .. math:: V_{rb}(\phi_{ijkl}) = \sum_{n=0}^5 C_n( \cos(\psi ))^n,

|  where :math:`\psi = \phi - 180^\circ`.
| **Note:** A conversion from one convention to another can be achieved
  by multiplying every coefficient :math:`\displaystyle C_n` by
  :math:`\displaystyle (-1)^n`.

An example of constants for :math:`C` is given in :numref:`Table %s <tab-crb>`.

.. _tab-crb:

.. table:: 
    Constants for Ryckaert-Bellemans potential (\ :math:`\mathrm{kJ mol}^{-1}`).
    :widths: auto
    :align: center

    +-------------+-------+-------------+--------+-------------+-------+
    | :math:`C_0` | 9.28  | :math:`C_2` | -13.12 | :math:`C_4` | 26.24 |
    +-------------+-------+-------------+--------+-------------+-------+
    | :math:`C_1` | 12.16 | :math:`C_3` | -3.06  | :math:`C_5` | -31.5 |
    +-------------+-------+-------------+--------+-------------+-------+


.. _fig-rbdih:

.. figure:: plots/f-rbs.*
   :width: 8.00000cm

   Ryckaert-Bellemans dihedral potential.

(**Note:** The use of this potential implies exclusion of LJ
interactions between the first and the last atom of the dihedral, and
:math:`\psi` is defined according to the “polymer convention”
(:math:`\psi_{trans}=0`).)

| The RB dihedral function can also be used to include Fourier dihedrals
  (see below):

  .. math::

     V_{rb} (\phi_{ijkl}) ~=~ \frac{1}{2} \left[F_1(1+\cos(\phi)) + F_2(
     1-\cos(2\phi)) + F_3(1+\cos(3\phi)) + F_4(1-\cos(4\phi))\right]

| Because of the equalities :math:`\cos(2\phi) = 2\cos^2(\phi) - 1`,
  :math:`\cos(3\phi) = 4\cos^3(\phi) - 3\cos(\phi)` and
  :math:`\cos(4\phi) = 8\cos^4(\phi) - 8\cos^2(\phi) + 1` one can
  translate the OPLS parameters to Ryckaert-Bellemans parameters as
  follows:

  .. math::

     \displaystyle
     \begin{array}{rcl}
     \displaystyle C_0&=&F_2 + \frac{1}{2} (F_1 + F_3)\\
     \displaystyle C_1&=&\frac{1}{2} (- F_1 + 3 \, F_3)\\
     \displaystyle C_2&=& -F_2 + 4 \, F_4\\
     \displaystyle C_3&=&-2 \, F_3\\
     \displaystyle C_4&=&-4 \, F_4\\
     \displaystyle C_5&=&0
     \end{array}

| with OPLS parameters in protein convention and RB parameters in
  polymer convention (this yields a minus sign for the odd powers of
  cos\ :math:`(\phi)`).
| **Note:** Mind the conversion from **kcal mol**\ :math:`^{-1}` for
  literature OPLS and RB parameters to **kJ mol**\ :math:`^{-1}` in
  |Gromacs|.

Proper dihedrals: Fourier function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The OPLS potential function is given as the first three
   :ref:`86 <refJorgensen1996>` or four \ :ref:`87 <refRobertson2015a>`
  cosine terms of a Fourier series. In |Gromacs| the four term function is
  implemented:

  .. math::

     V_{F} (\phi_{ijkl}) ~=~ \frac{1}{2} \left[C_1(1+\cos(\phi)) + C_2(
     1-\cos(2\phi)) + C_3(1+\cos(3\phi)) + C_4(1-\cos(4\phi))\right],

| Internally, |Gromacs| uses the Ryckaert-Bellemans code to compute
  Fourier dihedrals (see above), because this is more efficient.
| **Note:** Mind the conversion from *k*\ cal mol\ :math:`^{-1}` for
  literature OPLS parameters to **kJ mol**\ :math:`^{-1}` in |Gromacs|.

Proper dihedrals: Restricted torsion potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a manner very similar to the restricted bending potential (see
:ref:`ReB`), a restricted torsion/dihedral potential is introduced:

.. math:: V_{\rm ReT}(\phi_i) = \frac{1}{2} k_{\phi} \frac{(\cos\phi_i - \cos\phi_0)^2}{\sin^2\phi_i}
          :label: eqReT

with the advantages of being a function of :math:`\cos\phi` (no
problems taking the derivative of :math:`\sin\phi`) and of keeping the
torsion angle at only one minimum value. In this case, the factor
:math:`\sin^2\phi` does not allow the dihedral angle to move from the
[:math:`-180^{\circ}`:0] to [0::math:`180^{\circ}`] interval, i.e. it
cannot have maxima both at :math:`-\phi_0` and :math:`+\phi_0` maxima,
but only one of them. For this reason, all the dihedral angles of the
starting configuration should have their values in the desired angles
interval and the the equilibrium :math:`\phi_0` value should not be too
close to the interval limits (as for the restricted bending potential,
described in :ref:`ReB`, at least :math:`10^{\circ}` difference is
recommended).

Proper dihedrals: Combined bending-torsion potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the four particles forming the dihedral angle become collinear
(this situation will never happen in atomistic simulations, but it can
occur in coarse-grained simulations) the calculation of the torsion
angle and potential leads to numerical instabilities. One way to avoid
this is to use the restricted bending potential (see :ref:`ReB`) that
prevents the dihedral from reaching the :math:`180^{\circ}` value.

Another way is to disregard any effects of the dihedral becoming
ill-defined, keeping the dihedral force and potential calculation
continuous in entire angle range by coupling the torsion potential (in a
cosine form) with the bending potentials of the adjacent bending angles
in a unique expression:

.. math:: V_{\rm CBT}(\theta_{i-1}, \theta_i, \phi_i) = k_{\phi} \sin^3\theta_{i-1} \sin^3\theta_{i} \sum_{n=0}^4 { a_n \cos^n\phi_i}.
          :label: eqCBT

This combined bending-torsion (CBT) potential has been proposed
by \ :ref:`88 <refBulacuGiessen2005>` for polymer melt simulations and is
extensively described in \ :ref:`83 <refMonicaGoga2013>`.

This potential has two main advantages:

-  it does not only depend on the dihedral angle :math:`\phi_i` (between
   the :math:`i-2`, :math:`i-1`, :math:`i` and :math:`i+1` beads) but
   also on the bending angles :math:`\theta_{i-1}` and :math:`\theta_i`
   defined from three adjacent beads (:math:`i-2`, :math:`i-1` and
   :math:`i`, and :math:`i-1`, :math:`i` and :math:`i+1`, respectively).
   The two :math:`\sin^3\theta` pre-factors, tentatively suggested
   by \ :ref:`89 <refScottScheragator1966>` and theoretically discussed by
   \ :ref:`90 <refPaulingBond>`, cancel the torsion potential and force when either of the two
   bending angles approaches the value of :math:`180^\circ`.

-  its dependence on :math:`\phi_i` is expressed through a polynomial in
   :math:`\cos\phi_i` that avoids the singularities in
   :math:`\phi=0^\circ` or :math:`180^\circ` in calculating the
   torsional force.

These two properties make the CBT potential well-behaved for MD
simulations with weak constraints on the bending angles or even for
steered / non-equilibrium MD in which the bending and torsion angles
suffer major modifications. When using the CBT potential, the bending
potentials for the adjacent :math:`\theta_{i-1}` and :math:`\theta_i`
may have any form. It is also possible to leave out the two angle
bending terms (:math:`\theta_{i-1}` and :math:`\theta_{i}`) completely.
:numref:`Fig. %s <fig-CBT>` illustrates the difference between a torsion potential
with and without the :math:`\sin^{3}\theta` factors (blue and gray
curves, respectively).

.. _fig-CBT:

.. figure:: plots/fig-04.*
   :width: 10.00000cm

   Blue: surface plot of the combined bending-torsion potential
   (:eq:`%s <eqCBT>` with :math:`k = 10` kJ mol\ :math:`^{-1}`,
   :math:`a_0=2.41`, :math:`a_1=-2.95`, :math:`a_2=0.36`,
   :math:`a_3=1.33`) when, for simplicity, the bending angles behave the
   same (:math:`\theta_1=\theta_2=\theta`). Gray: the same torsion
   potential without the :math:`\sin^{3}\theta` terms
   (Ryckaert-Bellemans type). :math:`\phi` is the dihedral angle.

Additionally, the derivative of :math:`V_{CBT}` with respect to the
Cartesian variables is straightforward:

.. math:: \frac{\partial V_{\rm CBT}(\theta_{i-1},\theta_i,\phi_i)} {\partial \vec r_{l}} = \frac{\partial V_{\rm CBT}}{\partial \theta_{i-1}} \frac{\partial \theta_{i-1}}{\partial \vec r_{l}} +
          \frac{\partial V_{\rm CBT}}{\partial \theta_{i  }} \frac{\partial \theta_{i  }}{\partial \vec r_{l}} +
          \frac{\partial V_{\rm CBT}}{\partial \phi_{i    }} \frac{\partial \phi_{i    }}{\partial \vec r_{l}}
          :label: eqforcecbt

The CBT is based on a cosine form without multiplicity, so it can only
be symmetrical around :math:`0^{\circ}`. To obtain an asymmetrical
dihedral angle distribution (e.g. only one maximum in
[:math:`-180^{\circ}`::math:`180^{\circ}`] interval), a standard torsion
potential such as harmonic angle or periodic cosine potentials should be
used instead of a CBT potential. However, these two forms have the
inconveniences of the force derivation (:math:`1/\sin\phi`) and of the
alignment of beads (:math:`\theta_i` or
:math:`\theta_{i-1} = 0^{\circ}, 180^{\circ}`). Coupling such
non-\ :math:`\cos\phi` potentials with :math:`\sin^3\theta` factors does
not improve simulation stability since there are cases in which
:math:`\theta` and :math:`\phi` are simultaneously :math:`180^{\circ}`.
The integration at this step would be possible (due to the cancelling of
the torsion potential) but the next step would be singular
(:math:`\theta` is not :math:`180^{\circ}` and :math:`\phi` is very
close to :math:`180^{\circ}`).

Tabulated bonded interaction functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| For full flexibility, any functional shape can be used for bonds,
  angles and dihedrals through user-supplied tabulated functions. The
  functional shapes are:

  .. math::

     \begin{aligned}
     V_b(r_{ij})      &=& k \, f^b_n(r_{ij}) \\
     V_a({\theta_{ijk}})       &=& k \, f^a_n({\theta_{ijk}}) \\
     V_d(\phi_{ijkl}) &=& k \, f^d_n(\phi_{ijkl})\end{aligned}

| where :math:`k` is a force constant in units of energy and :math:`f`
  is a cubic spline function; for details see :ref:`cubicspline`. For
  each interaction, the force constant :math:`k` and the table number
  :math:`n` are specified in the topology. There are two different types
  of bonds, one that generates exclusions (type 8) and one that does not
  (type 9). For details see :numref:`Table %s <tab-topfile2>`. The table files are
  supplied to the :ref:`mdrun <gmx mdrun>` program. After the table file name an
  underscore, the letter “b” for bonds, “a” for angles or “d” for
  dihedrals and the table number must be appended. For example, a
  tabulated bond with :math:`n=0` can be read from the file
  table\_b0.xvg. Multiple tables can be supplied simply by adding files
  with different values of :math:`n`, and are applied to the appropriate
  bonds, as specified in the topology (:numref:`Table %s <tab-topfile2>`). The format
  for the table files is three fixed-format columns of any suitable
  width. These columns must contain :math:`x`, :math:`f(x)`,
  :math:`-f'(x)`, and the values of :math:`x` should be uniformly
  spaced. Requirements for entries in the topology are given
  in :numref:`Table %s <tab-topfile2>`. The setup of the tables is as follows:
| **bonds**: :math:`x` is the distance in nm. For distances beyond the
  table length, :ref:`mdrun <gmx mdrun>` will quit with an error message.
| **angles**: :math:`x` is the angle in degrees. The table should go
  from 0 up to and including 180 degrees; the derivative is taken in
  degrees.
| **dihedrals**: :math:`x` is the dihedral angle in degrees. The table
  should go from -180 up to and including 180 degrees; the IUPAC/IUB
  convention is used, *i.e.* zero is cis, the derivative is taken in
  degrees.

Restraints
----------

Special potentials are used for imposing restraints on the motion of the
system, either to avoid disastrous deviations, or to include knowledge
from experimental data. In either case they are not really part of the
force field and the reliability of the parameters is not important. The
potential forms, as implemented in |Gromacs|, are mentioned just for the
sake of completeness. Restraints and constraints refer to quite
different algorithms in |Gromacs|.

.. _positionrestraint:

Position restraints
~~~~~~~~~~~~~~~~~~~

These are used to restrain particles to fixed reference positions
:math:`\mathbf{R}_i`. They can be used during
equilibration in order to avoid drastic rearrangements of critical parts
(*e.g.* to restrain motion in a protein that is subjected to large
solvent forces when the solvent is not yet equilibrated). Another
application is the restraining of particles in a shell around a region
that is simulated in detail, while the shell is only approximated
because it lacks proper interaction from missing particles outside the
shell. Restraining will then maintain the integrity of the inner part.
For spherical shells, it is a wise procedure to make the force constant
depend on the radius, increasing from zero at the inner boundary to a
large value at the outer boundary. This feature has not, however, been
implemented in |Gromacs|.

The following form is used:

.. math:: V_{pr}(\mathbf{r}_i) = {\frac{1}{2}}k_{pr}|\mathbf{r}_i-\mathbf{R}_i|^2

The potential is plotted in :numref:`Fig. %s <fig-positionrestraint>`.

.. _fig-positionrestraint:

.. figure:: plots/f-pr.*
   :width: 8.00000cm

   Position restraint potential.

The potential form can be rewritten without loss of generality as:

.. math:: V_{pr}(\mathbf{r}_i) = {\frac{1}{2}} \left[ k_{pr}^x (x_i-X_i)^2 ~{\hat{\bf x}} + k_{pr}^y (y_i-Y_i)^2 ~{\hat{\bf y}} + k_{pr}^z (z_i-Z_i)^2 ~{\hat{\bf z}}\right]

Now the forces are:

.. math::

   \begin{array}{rcl}
   F_i^x &=& -k_{pr}^x~(x_i - X_i) \\
   F_i^y &=& -k_{pr}^y~(y_i - Y_i) \\
   F_i^z &=& -k_{pr}^z~(z_i - Z_i)
   \end{array}

Using three different force constants the position restraints can be
turned on or off in each spatial dimension; this means that atoms can be
harmonically restrained to a plane or a line. Position restraints are
applied to a special fixed list of atoms. Such a list is usually
generated by the :ref:`pdb2gmx <gmx pdb2gmx>` program.

Flat-bottomed position restraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Flat-bottomed position restraints can be used to restrain particles to
part of the simulation volume. No force acts on the restrained particle
within the flat-bottomed region of the potential, however a harmonic
force acts to move the particle to the flat-bottomed region if it is
outside it. It is possible to apply normal and flat-bottomed position
restraints on the same particle (however, only with the same reference
position :math:`\mathbf{R}_i`). The following general
potential is used (:numref:`Figure %s <fig-fbposres>` A):

.. math:: V_\mathrm{fb}(\mathbf{r}_i) = \frac{1}{2}k_\mathrm{fb} [d_g(\mathbf{r}_i;\mathbf{R}_i) - r_\mathrm{fb}]^2\,H[d_g(\mathbf{r}_i;\mathbf{R}_i) - r_\mathrm{fb}],

where :math:`\mathbf{R}_i` is the reference position,
:math:`r_\mathrm{fb}` is the distance from the center with a flat
potential, :math:`k_\mathrm{fb}` the force constant, and :math:`H` is
the Heaviside step function. The distance
:math:`d_g(\mathbf{r}_i;\mathbf{R}_i)` from
the reference position depends on the geometry :math:`g` of the
flat-bottomed potential.

.. _fig-fbposres:

.. figure:: plots/fbposres.*
   :width: 10.00000cm

   Flat-bottomed position restraint potential. (A) Not inverted, (B)
   inverted.

| The following geometries for the flat-bottomed potential are
  supported:

| **Sphere** (:math:`g =1`): The
  particle is kept in a sphere of given radius. The force acts towards
  the center of the sphere. The following distance calculation is used:

  .. math:: d_g(\mathbf{r}_i;\mathbf{R}_i) = | \mathbf{r}_i-\mathbf{R}_i |

| **Cylinder** (:math:`g=6,7,8`): The particle is kept in a cylinder of
  given radius parallel to the :math:`x` (:math:`g=6`), :math:`y`
  (:math:`g=7`), or :math:`z`-axis (:math:`g=8`). For backwards
  compatibility, setting :math:`g=2` is mapped to :math:`g=8` in the
  code so that old :ref:`tpr` files and topologies work. The
  force from the flat-bottomed potential acts towards the axis of the
  cylinder. The component of the force parallel to the cylinder axis is
  zero. For a cylinder aligned along the :math:`z`-axis:

  .. math:: d_g(\mathbf{r}_i;\mathbf{R}_i) = \sqrt{ (x_i-X_i)^2 + (y_i - Y_i)^2 }

| **Layer** (:math:`g=3,4,5`): The particle is kept in a layer defined
  by the thickness and the normal of the layer. The layer normal can be
  parallel to the :math:`x`, :math:`y`, or :math:`z`-axis. The force
  acts parallel to the layer normal.

  .. math::

     d_g(\mathbf{r}_i;\mathbf{R}_i) = |x_i-X_i|, \;\;\;\mbox{or}\;\;\; 
      d_g(\mathbf{r}_i;\mathbf{R}_i) = |y_i-Y_i|, \;\;\;\mbox{or}\;\;\; 
     d_g(\mathbf{r}_i;\mathbf{R}_i) = |z_i-Z_i|.

It is possible to apply multiple independent flat-bottomed position
restraints of different geometry on one particle. For example, applying
a cylinder and a layer in :math:`z` keeps a particle within a disk.
Applying three layers in :math:`x`, :math:`y`, and :math:`z` keeps the
particle within a cuboid.

In addition, it is possible to invert the restrained region with the
unrestrained region, leading to a potential that acts to keep the
particle *outside* of the volume defined by
:math:`\mathbf{R}_i`, :math:`g`, and
:math:`r_\mathrm{fb}`. That feature is switched on by defining a
negative :math:`r_\mathrm{fb}` in the topology. The following potential
is used (:numref:`Figure %s <fig-fbposres>` B):

.. math::

   V_\mathrm{fb}^{\mathrm{inv}}(\mathbf{r}_i) = \frac{1}{2}k_\mathrm{fb}
     [d_g(\mathbf{r}_i;\mathbf{R}_i) - | r_\mathrm{fb} | ]^2\,
     H[ -(d_g(\mathbf{r}_i;\mathbf{R}_i) - | r_\mathrm{fb} | )].

Angle restraints
~~~~~~~~~~~~~~~~

These are used to restrain the angle between two pairs of particles or
between one pair of particles and the :math:`z`-axis. The functional
form is similar to that of a proper dihedral. For two pairs of atoms:

.. math::

   V_{ar}(\mathbf{r}_i,\mathbf{r}_j,\mathbf{r}_k,\mathbf{r}_l)
           = k_{ar}(1 - \cos(n (\theta - \theta_0))
           )
   ,~~~~\mbox{where}~~
   \theta = \arccos\left(\frac{\mathbf{r}_j -\mathbf{r}_i}{\|\mathbf{r}_j -\mathbf{r}_i\|}
    \cdot \frac{\mathbf{r}_l -\mathbf{r}_k}{\|\mathbf{r}_l -\mathbf{r}_k\|} \right)

For one pair of atoms and the :math:`z`-axis:

.. math::

   V_{ar}(\mathbf{r}_i,\mathbf{r}_j) = k_{ar}(1 - \cos(n (\theta - \theta_0))
           )
   ,~~~~\mbox{where}~~
   \theta = \arccos\left(\frac{\mathbf{r}_j -\mathbf{r}_i}{\|\mathbf{r}_j -\mathbf{r}_i\|}
    \cdot \left( \begin{array}{c} 0 \\ 0 \\ 1 \\ \end{array} \right) \right)

A multiplicity (:math:`n`) of 2 is useful when you do not want to
distinguish between parallel and anti-parallel vectors. The equilibrium
angle :math:`\theta` should be between 0 and 180 degrees for
multiplicity 1 and between 0 and 90 degrees for multiplicity 2.

.. _dihedralrestraint:

Dihedral restraints
~~~~~~~~~~~~~~~~~~~

These are used to restrain the dihedral angle :math:`\phi` defined by
four particles as in an improper dihedral (sec. :ref:`imp`) but with a
slightly modified potential. Using:

.. math:: \phi' = \left(\phi-\phi_0\right) ~{\rm MOD}~ 2\pi
          :label: eqndphi

where :math:`\phi_0` is the reference angle, the potential is defined
as:

.. math:: V_{dihr}(\phi') ~=~ \left\{
          \begin{array}{lcllll}
          {\frac{1}{2}}k_{dihr}(\phi'-\phi_0-\Delta\phi)^2      
                          &\mbox{for}&     \phi' & >   & \Delta\phi       \\[1.5ex]
          0               &\mbox{for}&     \phi' & \le & \Delta\phi       \\[1.5ex]
          \end{array}\right.
          :label: eqndihre

where :math:`\Delta\phi` is a user defined angle and :math:`k_{dihr}`
is the force constant. **Note** that in the input in topology files,
angles are given in degrees and force constants in
kJ/mol/rad\ :math:`^2`.

.. _distancerestraint:

Distance restraints
~~~~~~~~~~~~~~~~~~~

Distance restraints add a penalty to the potential when the distance
between specified pairs of atoms exceeds a threshold value. They are
normally used to impose experimental restraints from, for instance,
experiments in nuclear magnetic resonance (NMR), on the motion of the
system. Thus, MD can be used for structure refinement using NMR data. In
|Gromacs| there are three ways to impose restraints on pairs of atoms:

-  Simple harmonic restraints: use ``[ bonds ]`` type 6 (see sec. :ref:`excl`).

-  Piecewise linear/harmonic restraints: ``[ bonds ]`` type
   10.

-  Complex NMR distance restraints, optionally with pair, time and/or
   ensemble averaging.

The last two options will be detailed now.

The potential form for distance restraints is quadratic below a
specified lower bound and between two specified upper bounds, and linear
beyond the largest bound (see :numref:`Fig. %s <fig-dist>`).

.. math:: V_{dr}(r_{ij}) ~=~ \left\{
          \begin{array}{lcllllll}
          {\frac{1}{2}}k_{dr}(r_{ij}-r_0)^2      
                          &\mbox{for}&     &     & r_{ij} & < & r_0       \\[1.5ex]
          0               &\mbox{for}& r_0 & \le & r_{ij} & < & r_1       \\[1.5ex]
          {\frac{1}{2}}k_{dr}(r_{ij}-r_1)^2      
                          &\mbox{for}& r_1 & \le & r_{ij} & < & r_2       \\[1.5ex]
          {\frac{1}{2}}k_{dr}(r_2-r_1)(2r_{ij}-r_2-r_1)  
                          &\mbox{for}& r_2 & \le & r_{ij} &   &
          \end{array}\right.
          :label: eqndisre

.. _fig-dist:

.. figure:: plots/f-dr.*
   :width: 8.00000cm

   Distance Restraint potential.

The forces are

.. math::

   \mathbf{F}_i~=~ \left\{
   \begin{array}{lcllllll}
   -k_{dr}(r_{ij}-r_0)\frac{\mathbf{r}_ij}{r_{ij}} 
                   &\mbox{for}&     &     & r_{ij} & < & r_0       \\[1.5ex]
   0               &\mbox{for}& r_0 & \le & r_{ij} & < & r_1       \\[1.5ex]
   -k_{dr}(r_{ij}-r_1)\frac{\mathbf{r}_ij}{r_{ij}} 
                   &\mbox{for}& r_1 & \le & r_{ij} & < & r_2       \\[1.5ex]
   -k_{dr}(r_2-r_1)\frac{\mathbf{r}_ij}{r_{ij}}    
                   &\mbox{for}& r_2 & \le & r_{ij} &   &
   \end{array} \right.

For restraints not derived from NMR data, this functionality will
usually suffice and a section of ``[ bonds ]`` type 10 can be used to apply individual
restraints between pairs of atoms, see :ref:`topfile`. For applying
restraints derived from NMR measurements, more complex functionality
might be required, which is provided through the ``[ distance_restraints ]`` section and is
described below.

Time averaging
^^^^^^^^^^^^^^

Distance restraints based on instantaneous distances can potentially
reduce the fluctuations in a molecule significantly. This problem can be
overcome by restraining to a *time averaged*
distance \ :ref:`91 <refTorda89>`. The forces with time averaging are:

.. math::

   \mathbf{F}_i~=~ \left\{
   \begin{array}{lcllllll}
   -k^a_{dr}(\bar{r}_{ij}-r_0)\frac{\mathbf{r}_ij}{r_{ij}}   
                   &\mbox{for}&     &     & \bar{r}_{ij} & < & r_0 \\[1.5ex]
   0               &\mbox{for}& r_0 & \le & \bar{r}_{ij} & < & r_1 \\[1.5ex]
   -k^a_{dr}(\bar{r}_{ij}-r_1)\frac{\mathbf{r}_ij}{r_{ij}}   
                   &\mbox{for}& r_1 & \le & \bar{r}_{ij} & < & r_2 \\[1.5ex]
   -k^a_{dr}(r_2-r_1)\frac{\mathbf{r}_ij}{r_{ij}}    
                   &\mbox{for}& r_2 & \le & \bar{r}_{ij} &   &
   \end{array} \right.

where :math:`\bar{r}_{ij}` is given by an exponential running average
with decay time :math:`\tau`:

.. math:: \bar{r}_{ij} ~=~ < r_{ij}^{-3} >^{-1/3}
          :label: eqnrav

The force constant :math:`k^a_{dr}` is switched on slowly to compensate
for the lack of history at the beginning of the simulation:

.. math:: k^a_{dr} = k_{dr} \left(1-\exp\left(-\frac{t}{\tau}\right)\right)

Because of the time averaging, we can no longer speak of a distance
restraint potential.

This way an atom can satisfy two incompatible distance restraints *on
average* by moving between two positions. An example would be an amino
acid side-chain that is rotating around its :math:`\chi` dihedral angle,
thereby coming close to various other groups. Such a mobile side chain
can give rise to multiple NOEs that can not be fulfilled by a single
structure.

The computation of the time averaged distance in the
:ref:`mdrun <gmx mdrun>` program is done in the following fashion:

.. math:: \begin{array}{rcl}
          \overline{r^{-3}}_{ij}(0)       &=& r_{ij}(0)^{-3}      \\
          \overline{r^{-3}}_{ij}(t)       &=& \overline{r^{-3}}_{ij}(t-\Delta t)~\exp{\left(-\frac{\Delta t}{\tau}\right)} + r_{ij}(t)^{-3}\left[1-\exp{\left(-\frac{\Delta t}{\tau}\right)}\right]
          \end{array}
          :label: eqnravdisre

When a pair is within the bounds, it can still feel a force because the
time averaged distance can still be beyond a bound. To prevent the
protons from being pulled too close together, a mixed approach can be
used. In this approach, the penalty is zero when the instantaneous
distance is within the bounds, otherwise the violation is the square
root of the product of the instantaneous violation and the time averaged
violation:

.. math::

   \mathbf{F}_i~=~ \left\{
   \begin{array}{lclll}
   k^a_{dr}\sqrt{(r_{ij}-r_0)(\bar{r}_{ij}-r_0)}\frac{\mathbf{r}_ij}{r_{ij}}   
       & \mbox{for} & r_{ij} < r_0 & \mbox{and} & \bar{r}_{ij} < r_0 \\[1.5ex]
   -k^a _{dr} \,
     \mbox{min}\left(\sqrt{(r_{ij}-r_1)(\bar{r}_{ij}-r_1)},r_2-r_1\right)
     \frac{\mathbf{r}_ij}{r_{ij}}   
       & \mbox{for} & r_{ij} > r_1 & \mbox{and} & \bar{r}_{ij} > r_1 \\[1.5ex]
   0               &\mbox{otherwise}
   \end{array} \right.

Averaging over multiple pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes it is unclear from experimental data which atom pair gives
rise to a single NOE, in other occasions it can be obvious that more
than one pair contributes due to the symmetry of the system, *e.g.* a
methyl group with three protons. For such a group, it is not possible to
distinguish between the protons, therefore they should all be taken into
account when calculating the distance between this methyl group and
another proton (or group of protons). Due to the physical nature of
magnetic resonance, the intensity of the NOE signal is inversely
proportional to the sixth power of the inter-atomic distance. Thus, when
combining atom pairs, a fixed list of :math:`N` restraints may be taken
together, where the apparent “distance” is given by:

.. math:: r_N(t) = \left [\sum_{n=1}^{N} \bar{r}_{n}(t)^{-6} \right]^{-1/6}
          :label: eqnrsix

where we use :math:`r_{ij}` or :eq:`eqn. %s <eqnrav>` for the
:math:`\bar{r}_{n}`. The :math:`r_N` of the instantaneous and
time-averaged distances can be combined to do a mixed restraining, as
indicated above. As more pairs of protons contribute to the same NOE
signal, the intensity will increase, and the summed “distance” will be
shorter than any of its components due to the reciprocal summation.

There are two options for distributing the forces over the atom pairs.
In the conservative option, the force is defined as the derivative of
the restraint potential with respect to the coordinates. This results in
a conservative potential when time averaging is not used. The force
distribution over the pairs is proportional to :math:`r^{-6}`. This
means that a close pair feels a much larger force than a distant pair,
which might lead to a molecule that is “too rigid.” The other option is
an equal force distribution. In this case each pair feels :math:`1/N` of
the derivative of the restraint potential with respect to :math:`r_N`.
The advantage of this method is that more conformations might be
sampled, but the non-conservative nature of the forces can lead to local
heating of the protons.

It is also possible to use *ensemble averaging* using multiple (protein)
molecules. In this case the bounds should be lowered as in:

.. math::

   \begin{array}{rcl}
   r_1     &~=~&   r_1 * M^{-1/6}  \\
   r_2     &~=~&   r_2 * M^{-1/6}
   \end{array}

where :math:`M` is the number of molecules. The |Gromacs| preprocessor
:ref:`grompp <gmx grompp>` can do this automatically when the appropriate
option is given. The resulting “distance” is then used to calculate the
scalar force according to:

.. math::

   \mathbf{F}_i~=~\left\{
   \begin{array}{rcl}
   ~& 0 \hspace{4cm}  & r_{N} < r_1         \\
    & k_{dr}(r_{N}-r_1)\frac{\mathbf{r}_ij}{r_{ij}} & r_1 \le r_{N} < r_2 \\
    & k_{dr}(r_2-r_1)\frac{\mathbf{r}_ij}{r_{ij}}    & r_{N} \ge r_2 
   \end{array} \right.

where :math:`i` and :math:`j` denote the atoms of all the pairs that
contribute to the NOE signal.

Using distance restraints
^^^^^^^^^^^^^^^^^^^^^^^^^

A list of distance restrains based on NOE data can be added to a
molecule definition in your topology file, like in the following
example:

::

    [ distance_restraints ]
    ; ai   aj   type   index   type'      low     up1     up2     fac
    10     16      1       0       1      0.0     0.3     0.4     1.0
    10     28      1       1       1      0.0     0.3     0.4     1.0
    10     46      1       1       1      0.0     0.3     0.4     1.0
    16     22      1       2       1      0.0     0.3     0.4     2.5
    16     34      1       3       1      0.0     0.5     0.6     1.0

In this example a number of features can be found. In columns ai and aj
you find the atom numbers of the particles to be restrained. The type
column should always be 1. As explained in  :ref:`distancerestraint`,
multiple distances can contribute to a single NOE signal. In the
topology this can be set using the index column. In our example, the
restraints 10-28 and 10-46 both have index 1, therefore they are treated
simultaneously. An extra requirement for treating restraints together is
that the restraints must be on successive lines, without any other
intervening restraint. The type’ column will usually be 1, but can be
set to 2 to obtain a distance restraint that will never be time- and
ensemble-averaged; this can be useful for restraining hydrogen bonds.
The columns ``low``, ``up1``, and
``up2`` hold the values of :math:`r_0`, :math:`r_1`, and
:math:`r_2` from  :eq:`eqn. %s <eqndisre>`. In some cases it
can be useful to have different force constants for some restraints;
this is controlled by the column ``fac``. The force constant
in the parameter file is multiplied by the value in the column
``fac`` for each restraint. Information for each restraint
is stored in the energy file and can be processed and plotted with
:ref:`gmx nmr`.

Orientation restraints
~~~~~~~~~~~~~~~~~~~~~~

This section describes how orientations between vectors, as measured in
certain NMR experiments, can be calculated and restrained in MD
simulations. The presented refinement methodology and a comparison of
results with and without time and ensemble averaging have been
published \ :ref:`92 <refHess2003>`.

Theory
^^^^^^

In an NMR experiment, orientations of vectors can be measured when a
molecule does not tumble completely isotropically in the solvent. Two
examples of such orientation measurements are residual dipolar couplings
(between two nuclei) or chemical shift anisotropies. An observable for a
vector :math:`\mathbf{r}_i` can be written as follows:

.. math:: \delta_i = \frac{2}{3} \mbox{tr}({{\mathbf S}}{{\mathbf D}}_i)

where :math:`{{\mathbf S}}` is the dimensionless order tensor of the
molecule. The tensor :math:`{{\mathbf D}}_i` is given by:

.. math:: {{\mathbf D}}_i = \frac{c_i}{\|\mathbf{r}_i\|^\alpha} \left(
          \begin{array}{lll}
          3 x x - 1 & 3 x y     & 3 x z     \\
          3 x y     & 3 y y - 1 & 3 y z     \\
          3 x z     & 3 y z     & 3 z z - 1 \\
          \end{array} \right)
          :label: eqnorientdef

.. math::

   \mbox{with:} \quad 
   x=\frac{r_{i,x}}{\|\mathbf{r}_i\|}, \quad
   y=\frac{r_{i,y}}{\|\mathbf{r}_i\|}, \quad 
   z=\frac{r_{i,z}}{\|\mathbf{r}_i\|}

For a dipolar coupling :math:`\mathbf{r}_i` is the vector
connecting the two nuclei, :math:`\alpha=3` and the constant :math:`c_i`
is given by:

.. math:: c_i = \frac{\mu_0}{4\pi} \gamma_1^i \gamma_2^i \frac{\hbar}{4\pi}

where :math:`\gamma_1^i` and :math:`\gamma_2^i` are the gyromagnetic
ratios of the two nuclei.

The order tensor is symmetric and has trace zero. Using a rotation
matrix :math:`{\mathbf T}` it can be transformed into the following
form:

.. math::

   {\mathbf T}^T {{\mathbf S}}{\mathbf T} = s \left( \begin{array}{ccc}
   -\frac{1}{2}(1-\eta) & 0                    & 0 \\
   0                    & -\frac{1}{2}(1+\eta) & 0 \\
   0                    & 0                    & 1
   \end{array} \right)

where :math:`-1 \leq s \leq 1` and :math:`0 \leq \eta \leq 1`.
:math:`s` is called the order parameter and :math:`\eta` the asymmetry
of the order tensor :math:`{{\mathbf S}}`. When the molecule tumbles
isotropically in the solvent, :math:`s` is zero, and no orientational
effects can be observed because all :math:`\delta_i` are zero.

Calculating orientations in a simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For reasons which are explained below, the :math:`{{\mathbf D}}`
matrices are calculated which respect to a reference orientation of the
molecule. The orientation is defined by a rotation matrix
:math:`{{\mathbf R}}`, which is needed to least-squares fit the current
coordinates of a selected set of atoms onto a reference conformation.
The reference conformation is the starting conformation of the
simulation. In case of ensemble averaging, which will be treated later,
the structure is taken from the first subsystem. The calculated
:math:`{{\mathbf D}}_i^c` matrix is given by:

.. math:: {{\mathbf D}}_i^c(t) = {{\mathbf R}}(t) {{\mathbf D}}_i(t) {{\mathbf R}}^T(t)
          :label: eqnDrot

The calculated orientation for vector :math:`i` is given by:

.. math:: \delta^c_i(t) = \frac{2}{3} \mbox{tr}({{\mathbf S}}(t){{\mathbf D}}_i^c(t))

The order tensor :math:`{{\mathbf S}}(t)` is usually unknown. A
reasonable choice for the order tensor is the tensor which minimizes the
(weighted) mean square difference between the calculated and the
observed orientations:

.. math:: MSD(t) = \left(\sum_{i=1}^N w_i\right)^{-1} \sum_{i=1}^N w_i (\delta_i^c (t) -\delta_i^{exp})^2
          :label: eqnSmsd

To properly combine different types of measurements, the unit of
:math:`w_i` should be such that all terms are dimensionless. This means
the unit of :math:`w_i` is the unit of :math:`\delta_i` to the power
:math:`-2`. **Note** that scaling all :math:`w_i` with a constant factor
does not influence the order tensor.

Time averaging
^^^^^^^^^^^^^^

Since the tensors :math:`{{\mathbf D}}_i` fluctuate rapidly in time,
much faster than can be observed in an experiment, they should be
averaged over time in the simulation. However, in a simulation the time
and the number of copies of a molecule are limited. Usually one can not
obtain a converged average of the :math:`{{\mathbf D}}_i` tensors over
all orientations of the molecule. If one assumes that the average
orientations of the :math:`\mathbf{r}_i` vectors within
the molecule converge much faster than the tumbling time of the
molecule, the tensor can be averaged in an axis system that rotates with
the molecule, as expressed by :eq:`equation %s <eqnDrot>`). The time-averaged
tensors are calculated using an exponentially decaying memory function:

.. math::

   {{\mathbf D}}^a_i(t) = \frac{\displaystyle
   \int_{u=t_0}^t {{\mathbf D}}^c_i(u) \exp\left(-\frac{t-u}{\tau}\right)\mbox{d} u
   }{\displaystyle
   \int_{u=t_0}^t \exp\left(-\frac{t-u}{\tau}\right)\mbox{d} u
   }

Assuming that the order tensor :math:`{{\mathbf S}}` fluctuates slower
than the :math:`{{\mathbf D}}_i`, the time-averaged orientation can be
calculated as:

.. math:: \delta_i^a(t) = \frac{2}{3} \mbox{tr}({{\mathbf S}}(t) {{\mathbf D}}_i^a(t))

where the order tensor :math:`{{\mathbf S}}(t)` is calculated using
expression :eq:`%s <eqnSmsd>` with :math:`\delta_i^c(t)` replaced by
:math:`\delta_i^a(t)`.

Restraining
^^^^^^^^^^^

The simulated structure can be restrained by applying a force
proportional to the difference between the calculated and the
experimental orientations. When no time averaging is applied, a proper
potential can be defined as:

.. math:: V = \frac{1}{2} k \sum_{i=1}^N w_i (\delta_i^c (t) -\delta_i^{exp})^2

where the unit of :math:`k` is the unit of energy. Thus the effective
force constant for restraint :math:`i` is :math:`k w_i`. The forces are
given by minus the gradient of :math:`V`. The force
:math:`\mathbf{F}\!_i` working on vector
:math:`\mathbf{r}_i` is:

.. math::

   \begin{aligned}
   \mathbf{F}\!_i(t) 
   & = & - \frac{\mbox{d} V}{\mbox{d}\mathbf{r}_i} \\
   & = & -k w_i (\delta_i^c (t) -\delta_i^{exp}) \frac{\mbox{d} \delta_i (t)}{\mbox{d}\mathbf{r}_i} \\
   & = & -k w_i (\delta_i^c (t) -\delta_i^{exp})
   \frac{2 c_i}{\|\mathbf{r}\|^{2+\alpha}} \left(2 {{\mathbf R}}^T {{\mathbf S}}{{\mathbf R}}\mathbf{r}_i - \frac{2+\alpha}{\|\mathbf{r}\|^2} \mbox{tr}({{\mathbf R}}^T {{\mathbf S}}{{\mathbf R}}\mathbf{r}_i \mathbf{r}_i^T) \mathbf{r}_i \right)\end{aligned}

Ensemble averaging
^^^^^^^^^^^^^^^^^^

Ensemble averaging can be applied by simulating a system of :math:`M`
subsystems that each contain an identical set of orientation restraints.
The systems only interact via the orientation restraint potential which
is defined as:

.. math::

   V = M \frac{1}{2} k \sum_{i=1}^N w_i 
   \langle \delta_i^c (t) -\delta_i^{exp} \rangle^2

The force on vector :math:`\mathbf{r}_{i,m}` in subsystem
:math:`m` is given by:

.. math::

   \mathbf{F}\!_{i,m}(t) = - \frac{\mbox{d} V}{\mbox{d}\mathbf{r}_{i,m}} =
   -k w_i \langle \delta_i^c (t) -\delta_i^{exp} \rangle \frac{\mbox{d} \delta_{i,m}^c (t)}{\mbox{d}\mathbf{r}_{i,m}} \\

Time averaging
^^^^^^^^^^^^^^

When using time averaging it is not possible to define a potential. We
can still define a quantity that gives a rough idea of the energy stored
in the restraints:

.. math::

   V = M \frac{1}{2} k^a \sum_{i=1}^N w_i 
   \langle \delta_i^a (t) -\delta_i^{exp} \rangle^2

The force constant :math:`k_a` is switched on slowly to compensate for
the lack of history at times close to :math:`t_0`. It is exactly
proportional to the amount of average that has been accumulated:

.. math::

   k^a =
    k \, \frac{1}{\tau}\int_{u=t_0}^t \exp\left(-\frac{t-u}{\tau}\right)\mbox{d} u

What really matters is the definition of the force. It is chosen to be
proportional to the square root of the product of the time-averaged and
the instantaneous deviation. Using only the time-averaged deviation
induces large oscillations. The force is given by:

.. math::

   \mathbf{F}\!_{i,m}(t) =
   \left\{ \begin{array}{ll}
   0 & \quad \mbox{for} \quad a\, b \leq 0 \\
   \displaystyle
   k^a w_i \frac{a}{|a|} \sqrt{a\, b} \, \frac{\mbox{d} \delta_{i,m}^c (t)}{\mbox{d}\mathbf{r}_{i,m}}
   & \quad \mbox{for} \quad a\, b > 0 
   \end{array}
   \right.

.. math::

   \begin{aligned}
   a &=& \langle \delta_i^a (t) -\delta_i^{exp} \rangle \\
   b &=& \langle \delta_i^c (t) -\delta_i^{exp} \rangle\end{aligned}

Using orientation restraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Orientation restraints can be added to a molecule definition in the
topology file in the section ``[ orientation_restraints ]``.
Here we give an example section containing five N-H residual dipolar
coupling restraints:

::

    [ orientation_restraints ]
    ; ai   aj  type  exp.  label  alpha    const.     obs.   weight
    ;                                Hz      nm^3       Hz    Hz^-2
      31   32     1     1      3      3     6.083    -6.73      1.0
      43   44     1     1      4      3     6.083    -7.87      1.0
      55   56     1     1      5      3     6.083    -7.13      1.0
      65   66     1     1      6      3     6.083    -2.57      1.0
      73   74     1     1      7      3     6.083    -2.10      1.0

The unit of the observable is Hz, but one can choose any other unit. In
columns ``ai`` and ``aj`` you find the atom numbers of the particles to be
restrained. The ``type`` column should always be 1. The ``exp.`` column denotes
the experiment number, starting at 1. For each experiment a separate
order tensor :math:`{{\mathbf S}}` is optimized. The label should be a
unique number larger than zero for each restraint. The ``alpha`` column
contains the power :math:`\alpha` that is used in
:eq:`equation %s <eqnorientdef>`) to calculate the orientation. The ``const.`` column
contains the constant :math:`c_i` used in the same equation. The
constant should have the unit of the observable times
nm\ :math:`^\alpha`. The column ``obs.`` contains the observable, in any
unit you like. The last column contains the weights :math:`w_i`; the
unit should be the inverse of the square of the unit of the observable.

Some parameters for orientation restraints can be specified in the
:ref:`grompp <gmx grompp>` :ref:`mdp` file, for a study of the effect of different
force constants and averaging times and ensemble averaging see \ :ref:`92 <refHess2003>`.
Information for each restraint is stored in the energy
file and can be processed and plotted with :ref:`gmx nmr`.

Polarization
------------

Polarization can be treated by |Gromacs| by attaching shell (Drude)
particles to atoms and/or virtual sites. The energy of the shell
particle is then minimized at each time step in order to remain on the
Born-Oppenheimer surface.

Simple polarization
~~~~~~~~~~~~~~~~~~~

This is implemented as a harmonic potential with equilibrium distance 0.
The input given in the topology file is the polarizability
:math:`\alpha` (in |Gromacs| units) as follows:

::

    [ polarization ]
    ; Atom i  j  type  alpha
    1         2  1     0.001

in this case the polarizability volume is 0.001 nm\ :math:`^3` (or 1
Å\ :math:`^3`). In order to compute the harmonic force constant
:math:`k_{cs}` (where :math:`cs` stands for core-shell), the following
is used \ :ref:`45 <refMaaren2001a>`:

.. math:: k_{cs} ~=~ \frac{q_s^2}{\alpha}

where :math:`q_s` is the charge on the shell particle.

Anharmonic polarization
~~~~~~~~~~~~~~~~~~~~~~~

For the development of the Drude force field by Roux and
McKerell \ :ref:`93 <refLopes2013a>` it was found that some particles can
overpolarize and this was fixed by introducing a higher order term in
the polarization energy:

.. math::

   \begin{aligned}
   V_{pol} ~=& \frac{k_{cs}}{2} r_{cs}^2 & r_{cs} \le \delta \\
               =& \frac{k_{cs}}{2} r_{cs}^2 + k_{hyp} (r_{cs}-\delta)^4 & r_{cs} > \delta\end{aligned}

where :math:`\delta` is a user-defined constant that is set to 0.02 nm
for anions in the Drude force field \ :ref:`94 <refHYu2010>`. Since this
original introduction it has also been used in other atom
types \ :ref:`93 <refLopes2013a>`.

::

    [ polarization ]
    ;Atom i j    type   alpha (nm^3)    delta  khyp
    1       2       2       0.001786     0.02  16.736e8

The above force constant :math:`k_{hyp}` corresponds to
4\ :math:`\cdot`\ 10\ :math:`^8` kcal/mol/nm\ :math:`^4`, hence the
strange number.

Water polarization
~~~~~~~~~~~~~~~~~~

A special potential for water that allows anisotropic polarization of a
single shell particle \ :ref:`45 <refMaaren2001a>`.

Thole polarization
~~~~~~~~~~~~~~~~~~

Based on early work by Thole :ref:`95 <refThole81>`, Roux and coworkers
have implemented potentials for molecules like
ethanol \ :ref:`96 <refLamoureux2003a>`\ :ref:`98 <refNoskov2005a>`.
Within such molecules, there are intra-molecular interactions between
shell particles, however these must be screened because full Coulomb
would be too strong. The potential between two shell particles :math:`i`
and :math:`j` is:

.. math:: V_{thole} ~=~ \frac{q_i q_j}{r_{ij}}\left[1-\left(1+\frac{{\bar{r}_{ij}}}{2}\right){\rm exp}^{-{\bar{r}_{ij}}}\right]

**Note** that there is a sign error in Equation 1 of Noskov *et
al.* :ref:`98 <refNoskov2005a>`:

.. math:: {\bar{r}_{ij}}~=~ a\frac{r_{ij}}{(\alpha_i \alpha_j)^{1/6}}

where :math:`a` is a magic (dimensionless) constant, usually chosen to
be 2.6 \ :ref:`98 <refNoskov2005a>`; :math:`\alpha_i` and
:math:`\alpha_j` are the polarizabilities of the respective shell
particles.

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

.. math::

   \begin{aligned}
   V_b     &=&{\frac{1}{2}}\left[{(1-{\lambda})}k_b^A + 
                   {\lambda}k_b^B\right] \left[b - {(1-{\lambda})}b_0^A - {\lambda}b_0^B\right]^2  \\
   {\frac{\partial V_b}{\partial {\lambda}}}&=&{\frac{1}{2}}(k_b^B-k_b^A)
                   \left[b - {(1-{\lambda})}b_0^A + {\lambda}b_0^B\right]^2 + 
   		\nonumber\\
           & & \phantom{{\frac{1}{2}}}(b_0^A-b_0^B) \left[b - {(1-{\lambda})}b_0^A -{\lambda}b_0^B\right]
   		\left[{(1-{\lambda})}k_b^A + {\lambda}k_b^B \right]\end{aligned}

GROMOS-96 bonds and angles
~~~~~~~~~~~~~~~~~~~~~~~~~~

Fourth-power bond stretching and cosine-based angle potentials are
interpolated by linear interpolation of the force constant and the
equilibrium position. Formulas are not given here.

Proper dihedrals
~~~~~~~~~~~~~~~~

For the proper dihedrals, the equations are somewhat more complicated:

.. math::

   \begin{aligned}
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

**Note:** that the multiplicity :math:`n_{\phi}` can not be
parameterized because the function should remain periodic on the
interval :math:`[0,2\pi]`.

Tabulated bonded interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For tabulated bonded interactions only the force constant can
interpolated:

.. math::

   \begin{aligned}
         V  &=& ({(1-{\lambda})}k^A + {\lambda}k^B) \, f \\
   {\frac{\partial V}{\partial {\lambda}}} &=& (k^B - k^A) \, f\end{aligned}

Coulomb interaction
~~~~~~~~~~~~~~~~~~~

The Coulomb interaction between two particles of which the charge varies
with :math:`{\lambda}` is:

.. math::

   \begin{aligned}
   V_c &=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[{(1-{\lambda})}q_i^A q_j^A + {\lambda}\, q_i^B q_j^B\right] \\
   {\frac{\partial V_c}{\partial {\lambda}}}&=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[- q_i^A q_j^A + q_i^B q_j^B\right]\end{aligned}

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
m :math:`\mathbf{v}`, since that would result in the
sign of :math:`{\frac{\partial E_k}{\partial {\lambda}}}` being
incorrect \ :ref:`99 <refGunsteren98a>`):

.. math::

   \begin{aligned}
   E_k      &=&     {\frac{1}{2}}\frac{\mathbf{p}^2}{{(1-{\lambda})}m^A + {\lambda}m^B}        \\
   {\frac{\partial E_k}{\partial {\lambda}}}&=&    -{\frac{1}{2}}\frac{\mathbf{p}^2(m^B-m^A)}{({(1-{\lambda})}m^A + {\lambda}m^B)^2}\end{aligned}

after taking the derivative, we *can* insert
:math:`\mathbf{p}` =
m :math:`\mathbf{v}`, such that:

.. math:: {\frac{\partial E_k}{\partial {\lambda}}}~=~    -{\frac{1}{2}}\mathbf{v}^2(m^B-m^A)

Constraints
~~~~~~~~~~~

The constraints are formally part of the Hamiltonian, and therefore they
give a contribution to the free energy. In |Gromacs| this can be
calculated using the LINCS or the SHAKE algorithm. If we have
:math:`k = 1 \ldots K` constraint equations :math:`g_k` for LINCS, then

.. math:: g_k     =       | \mathbf{r}_{k} | - d_{k}

where :math:`\mathbf{r}_k` is the displacement vector
between two particles and :math:`d_k` is the constraint distance between
the two particles. We can express the fact that the constraint distance
has a :math:`{\lambda}` dependency by

.. math:: d_k     =       {(1-{\lambda})}d_{k}^A + {\lambda}d_k^B

Thus the :math:`{\lambda}`-dependent constraint equation is

.. math:: g_k     =       | \mathbf{r}_{k} | - \left({(1-{\lambda})}d_{k}^A + {\lambda}d_k^B\right).

The (zero) contribution :math:`G` to the Hamiltonian from the
constraints (using Lagrange multipliers :math:`\lambda_k`, which are
logically distinct from the free-energy :math:`{\lambda}`) is

.. math::

   \begin{aligned}
   G           &=&     \sum^K_k \lambda_k g_k    \\
   {\frac{\partial G}{\partial {\lambda}}}    &=&     \frac{\partial G}{\partial d_k} {\frac{\partial d_k}{\partial {\lambda}}} \\
               &=&     - \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}

For SHAKE, the constraint equations are

.. math:: g_k     =       \mathbf{r}_{k}^2 - d_{k}^2

with :math:`d_k` as before, so

.. math::

   \begin{aligned}
   {\frac{\partial G}{\partial {\lambda}}}    &=&     -2 \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}

Soft-core interactions
~~~~~~~~~~~~~~~~~~~~~~

.. _fig-softcore:

.. figure:: plots/softcore.*
   :height: 6.00000cm

   Soft-core interactions at :math:`{\lambda}=0.5`, with :math:`p=2` and
   :math:`C_6^A=C_{12}^A=C_6^B=C_{12}^B=1`.

In a free-energy calculation where particles grow out of nothing, or
particles disappear, using the the simple linear interpolation of the
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
its derivatives at :math:`r=0` is never reached:

.. math::

   \begin{aligned}
   V_{sc}(r) &=& {(1-{\lambda})}V^A(r_A) + {\lambda}V^B(r_B)
       \\
   r_A &=& \left(\alpha \sigma_A^6 {\lambda}^p + r^6 \right)^\frac{1}{6}
       \\
   r_B &=& \left(\alpha \sigma_B^6 {(1-{\lambda})}^p + r^6 \right)^\frac{1}{6}\end{aligned}

where :math:`V^A` and :math:`V^B` are the normal “hard core” Van der
Waals or electrostatic potentials in state A (:math:`{\lambda}=0`) and
state B (:math:`{\lambda}=1`) respectively, :math:`\alpha` is the
soft-core parameter (set with ``sc_alpha`` in the
:ref:`mdp` file), :math:`p` is the soft-core :math:`{\lambda}`
power (set with ``sc_power``), :math:`\sigma` is the radius
of the interaction, which is :math:`(C_{12}/C_6)^{1/6}` or an input
parameter (``sc_sigma``) when :math:`C_6` or :math:`C_{12}`
is zero.

For intermediate :math:`{\lambda}`, :math:`r_A` and :math:`r_B` alter
the interactions very little for :math:`r > \alpha^{1/6} \sigma` and
quickly switch the soft-core interaction to an almost constant value for
smaller :math:`r` (:numref:`Fig. %s <fig-softcore>`). The force is:

.. math::

   F_{sc}(r) = -\frac{\partial V_{sc}(r)}{\partial r} =
    {(1-{\lambda})}F^A(r_A) \left(\frac{r}{r_A}\right)^5 +
   {\lambda}F^B(r_B) \left(\frac{r}{r_B}\right)^5

where :math:`F^A` and :math:`F^B` are the “hard core” forces. The
contribution to the derivative of the free energy is:

.. math::

   \begin{aligned}
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

Recently, a new formulation of the soft-core approach has been derived
that in most cases gives lower and more even statistical variance than
the standard soft-core path described above \ :ref:`101 <refPham2011>`,
:ref:`102 <refPham2012>`. Specifically, we have:

.. math::

   \begin{aligned}
   V_{sc}(r) &=& {(1-{\lambda})}V^A(r_A) + {\lambda}V^B(r_B)
       \\
   r_A &=& \left(\alpha \sigma_A^{48} {\lambda}^p + r^{48} \right)^\frac{1}{48}
       \\
   r_B &=& \left(\alpha \sigma_B^{48} {(1-{\lambda})}^p + r^{48} \right)^\frac{1}{48}\end{aligned}

This “1-1-48” path is also implemented in |Gromacs|. Note that for this
path the soft core :math:`\alpha` should satisfy
:math:`0.001 < \alpha < 0.003`, rather than :math:`\alpha \approx
0.5`.

Methods
-------

Exclusions and 1-4 Interactions.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Atoms within a molecule that are close by in the chain, *i.e.* atoms
that are covalently bonded, or linked by one or two atoms are called
*first neighbors, second neighbors* and *third neighbors*, respectively
(see :numref:`Fig. %s <fig-chain>`). Since the interactions of atom **i** with atoms
**i+1** and **i+2** are mainly quantum mechanical, they can not be
modeled by a Lennard-Jones potential. Instead it is assumed that these
interactions are adequately modeled by a harmonic bond term or
constraint (**i, i+1**) and a harmonic angle term (**i, i+2**). The
first and second neighbors (atoms **i+1** and **i+2**) are therefore
*excluded* from the Lennard-Jones interaction list of atom **i**; atoms
**i+1** and **i+2** are called *exclusions* of atom **i**.

.. _fig-chain:

.. figure:: plots/chain.*
   :width: 8.00000cm

   Atoms along an alkane chain.

For third neighbors, the normal Lennard-Jones repulsion is sometimes
still too strong, which means that when applied to a molecule, the
molecule would deform or break due to the internal strain. This is
especially the case for carbon-carbon interactions in a
*cis*-conformation (*e.g.* *cis*-butane). Therefore, for some of these
interactions, the Lennard-Jones repulsion has been reduced in the GROMOS
force field, which is implemented by keeping a separate list of 1-4 and
normal Lennard-Jones parameters. In other force fields, such as
OPLS \ :ref:`103 <refJorgensen88>`, the standard Lennard-Jones
parameters are reduced by a factor of two, but in that case also the
dispersion (r:math:`^{-6}`) and the Coulomb interaction are scaled.
|Gromacs| can use either of these methods.

Charge Groups
~~~~~~~~~~~~~

In principle, the force calculation in MD is an :math:`O(N^2)` problem.
Therefore, we apply a cut-off for non-bonded force (NBF) calculations;
only the particles within a certain distance of each other are
interacting. This reduces the cost to :math:`O(N)` (typically
:math:`100N` to :math:`200N`) of the NBF. It also introduces an error,
which is, in most cases, acceptable, except when applying the cut-off
implies the creation of charges, in which case you should consider using
the lattice sum methods provided by |Gromacs|.

Consider a water molecule interacting with another atom. If we would
apply a plain cut-off on an atom-atom basis we might include the
atom-oxygen interaction (with a charge of :math:`-0.82`) without the
compensating charge of the protons, and as a result, induce a large
dipole moment over the system. Therefore, we have to keep groups of
atoms with total charge 0 together. These groups are called *charge
groups*. Note that with a proper treatment of long-range electrostatics
(e.g. particle-mesh Ewald (sec. :ref:`pme`), keeping charge groups
together is not required.

Treatment of Cut-offs in the group scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|Gromacs| is quite flexible in treating cut-offs, which implies there can
be quite a number of parameters to set. These parameters are set in the
input file for grompp. There are two sort of parameters that affect the
cut-off interactions; you can select which type of interaction to use in
each case, and which cut-offs should be used in the neighbor searching.

For both Coulomb and van der Waals interactions there are interaction
type selectors (termed vdwtype and coulombtype) and two parameters, for
a total of six non-bonded interaction parameters. See the User Guide for
a complete description of these parameters.

In the group cut-off scheme, all of the interaction functions in
:numref:`Table %s <tab-funcparm>` require that neighbor searching be done with a
radius at least as large as the :math:`r_c` specified for the functional
form, because of the use of charge groups. The extra radius is typically
of the order of 0.25 nm (roughly the largest distance between two atoms
in a charge group plus the distance a charge group can diffuse within
neighbor list updates).

.. |CPCOP| replace:: :math:`r_c`, :math:`{\varepsilon}_{r}`
.. |CRFP|  replace:: :math:`r_c`, :math:`{\varepsilon}_{rf}`
.. |CSHFP| replace:: :math:`r_1`, :math:`r_c`, :math:`{\varepsilon}_{r}`
.. |CSWFP| replace:: :math:`r_1`, :math:`r_c`, :math:`{\varepsilon}_{r}`
.. |VPCOP| replace:: :math:`r_c`
.. |VSHFP| replace:: :math:`r_1`, :math:`r_c`
.. |VSWFP| replace:: :math:`r_1`, :math:`r_c`

.. _tab-funcparm:

.. table:: Parameters for the different functional forms of the
           non-bonded interactions.

           +----------------------------+------------+
           | Type                       | Parameters |
           +=========+==================+============+
           | Coulomb | Plain cut-off    | |CPCOP|    |
           |         +------------------+------------+
           |         | Reaction field   | |CRFP|     |
           |         +------------------+------------+
           |         | Shift function   | |CSHFP|    |
           |         +------------------+------------+ 
           |         | Switch function  | |CSWFP|    | 
           +---------+------------------+------------+
           | VdW     | Plain cut-off    | |VPCOP|    |
           |         +------------------+------------+ 
           |         | Shift function   | |VSHFP|    |
           |         +------------------+------------+ 
           |         | Switch function  | |VSWFP|    | 
           +---------+------------------+------------+

.. _virtualsites:

Virtual interaction sites
-------------------------

Virtual interaction sites (called dummy atoms in
|Gromacs| versions before 3.3) can be used in |Gromacs| in a number of ways.
We write the position of the virtual site :math:`\mathbf{r}_s` as a function
of the positions of other particles
:math:`\mathbf{r}`\ :math:`_i`: :math:`\mathbf{r}_s =
f(\mathbf{r}_1..\mathbf{r}_n)`. The virtual site, which may carry charge or be
involved in other interactions, can now be used in the force
calculation. The force acting on the virtual site must be redistributed
over the particles with mass in a consistent way. A good way to do this
can be found in ref. \ :ref:`104 <refBerendsen84b>`. We can write the
potential energy as:

.. math:: V = V(\mathbf{r}_s,\mathbf{r}_1,\ldots,\mathbf{r}_n) = V^*(\mathbf{r}_1,\ldots,\mathbf{r}_n)

The force on the particle :math:`i` is then:

.. math::

   \mathbf{F}_i = -\frac{\partial V^*}{\partial \mathbf{r}_i} 
            = -\frac{\partial V}{\partial \mathbf{r}_i} - 
               \frac{\partial V}{\partial \mathbf{r}_s} 
               \frac{\partial \mathbf{r}_s}{\partial \mathbf{r}_i}
            = \mathbf{F}_i^{direct} + \mathbf{F}_i

The first term is the normal force. The second term is the force on
particle :math:`i` due to the virtual site, which can be written in
tensor notation:

.. math::  \mathbf{F}_i = \left[\begin{array}{ccc}
           {\displaystyle\frac{\partial x_s}{\partial x_i}} & {\displaystyle\frac{\partial y_s}{\partial x_i}} & {\displaystyle\frac{\partial z_s}{\partial x_i}} \\[1ex]
           {\displaystyle\frac{\partial x_s}{\partial y_i}} & {\displaystyle\frac{\partial y_s}{\partial y_i}} & {\displaystyle\frac{\partial z_s}{\partial y_i}} \\[1ex]
           {\displaystyle\frac{\partial x_s}{\partial z_i}} & {\displaystyle\frac{\partial y_s}{\partial z_i}} & {\displaystyle\frac{\partial z_s}{\partial z_i}} \end{array}\right]\mathbf{F}_{s}
           :label: eqnfvsite

where :math:`\mathbf{F}_{s}` is the force on the virtual site and
:math:`x_s`, :math:`y_s` and :math:`z_s` are the coordinates of the
virtual site. In this way, the total force and the total torque are
conserved \ :ref:`104 <refBerendsen84b>`.

The computation of the virial (:eq:`eqn. %s <eqnXi>`) is non-trivial when
virtual sites are used. Since the virial involves a summation over all
the atoms (rather than virtual sites), the forces must be redistributed
from the virtual sites to the atoms (using  :eq:`eqn. %s <eqnfvsite>`) *before*
computation of the virial. In some special cases where the forces on the
atoms can be written as a linear combination of the forces on the
virtual sites (types 2 and 3 below) there is no difference between
computing the virial before and after the redistribution of forces.
However, in the general case redistribution should be done first.

.. _fig-vsites:

.. figure:: plots/dummies.*
   :width: 15.00000cm

   The six different types of virtual site construction in . The
   constructing atoms are shown as black circles, the virtual sites in
   gray.

There are six ways to construct virtual sites from surrounding atoms in
|Gromacs|, which we classify by the number of constructing atoms. **Note**
that all site types mentioned can be constructed from types 3fd
(normalized, in-plane) and 3out (non-normalized, out of plane). However,
the amount of computation involved increases sharply along this list, so
we strongly recommended using the first adequate virtual site type that
will be sufficient for a certain purpose. :numref:`Fig. %s <fig-vsites>` depicts 6 of
the available virtual site constructions. The conceptually simplest
construction types are linear combinations:

.. math:: \mathbf{r}_s = \sum_{i=1}^N w_i \, \mathbf{r}_i

The force is then redistributed using the same weights:

.. math:: \mathbf{F}_i = w_i \, \mathbf{F}_{s}

The types of virtual sites supported in |Gromacs| are given in the list
below. Constructing atoms in virtual sites can be virtual sites
themselves, but only if they are higher in the list, i.e. virtual sites
can be constructed from “particles” that are simpler virtual sites.

-  As a linear combination of two atoms
   (:numref:`Fig. %s <fig-vsites>` 2):

   .. math:: w_i = 1 - a ~,~~ w_j = a

-  In this case the virtual site is on the line through atoms :math:`i`
   and :math:`j`.

-  As a linear combination of three atoms
   (:numref:`Fig. %s <fig-vsites>` 3):

   .. math:: w_i = 1 - a - b ~,~~ w_j = a ~,~~ w_k = b

-  In this case the virtual site is in the plane of the other three
   particles.

-  In the plane of three atoms, with a fixed distance
   (:numref:`Fig. %s <fig-vsites>` 3fd):

   .. math::

      \mathbf{r}_s ~=~ \mathbf{r}_i + b \frac{  \mathbf{r}_ij + a \mathbf{r}_{jk}  }
                                          { | \mathbf{r}_ij + a \mathbf{r}_{jk} | }

-  In this case the virtual site is in the plane of the other three
   particles at a distance of :math:`|b|` from :math:`i`. The force on
   particles :math:`i`, :math:`j` and :math:`k` due to the force on the
   virtual site can be computed as:

   .. math::

      \begin{array}{lcr}
              \mathbf{F}_i &=& \displaystyle \mathbf{F}_{s} - \gamma ( \mathbf{F}_is - \mathbf{p} ) \\[1ex]
              \mathbf{F}_j &=& \displaystyle (1-a)\gamma (\mathbf{F}_{s} - \mathbf{p})      \\[1ex]
              \mathbf{F}_k &=& \displaystyle a \gamma (\mathbf{F}_{s} - \mathbf{p})         \\
              \end{array}
              ~\mbox{~ where~ }~
              \begin{array}{c}
      \displaystyle \gamma = \frac{b}{ | \mathbf{r}_ij + a \mathbf{r}_{jk} | } \\[2ex]
      \displaystyle \mathbf{p} = \frac{ \mathbf{r}_{is} \cdot \mathbf{F}_{s} }
                            { \mathbf{r}_{is} \cdot \mathbf{r}_is } \mathbf{r}_is
              \end{array}

-  In the plane of three atoms, with a fixed angle and
   distance (:numref:`Fig. %s <fig-vsites>` 3fad):

   .. math:: \mathbf{r}_s ~=~ \mathbf{r}_i +
             d \cos \theta \frac{\mathbf{r}_ij}{ | \mathbf{r}_{ij} | } +
             d \sin \theta \frac{\mathbf{r}_\perp}{ | \mathbf{r}_\perp | }
             ~\mbox{~ where~ }~
             \mathbf{r}_\perp ~=~ \mathbf{r}_{jk} - 
             \frac{ \mathbf{r}_ij \cdot \mathbf{r}_{jk} }
             { \mathbf{r}_ij \cdot \mathbf{r}_{ij} }
             \mathbf{r}_ij
             :label: eqnvsite2fadF

-  In this case the virtual site is in the plane of the other three
   particles at a distance of :math:`|d|` from :math:`i` at an angle of
   :math:`\alpha` with :math:`\mathbf{r}_ij`. Atom
   :math:`k` defines the plane and the direction of the angle. **Note**
   that in this case :math:`b` and :math:`\alpha` must be specified,
   instead of :math:`a` and :math:`b` (see also sec. :ref:`vsitetop`).
   The force on particles :math:`i`, :math:`j` and :math:`k` due to the
   force on the virtual site can be computed as (with
   :math:`\mathbf{r}_\perp` as defined in
   :eq:`eqn. %s <eqnvsite2fadF>`):

   .. math::

      \begin{array}{c}
              \begin{array}{lclllll}
              \mathbf{F}_i &=& \mathbf{F}_{s} &-& 
                      \dfrac{d \cos \theta}{ | \mathbf{r}_ij | } \mathbf{F}_1 &+&
                      \dfrac{d \sin \theta}{ | \mathbf{r}_\perp | } \left( 
                      \dfrac{ \mathbf{r}_ij \cdot \mathbf{r}_{jk} }
                           { \mathbf{r}_ij \cdot \mathbf{r}_{ij} } \mathbf{F}_2     +
                      \mathbf{F}_3 \right)                                \\[3ex]
              \mathbf{F}_j &=& &&
                      \dfrac{d \cos \theta}{ | \mathbf{r}_ij | } \mathbf{F}_1 &-&
                      \dfrac{d \sin \theta}{ | \mathbf{r}_\perp | } \left(
                       \mathbf{F}_2 + 
                       \dfrac{ \mathbf{r}_ij \cdot \mathbf{r}_{jk} }
                              { \mathbf{r}_ij \cdot \mathbf{r}_{ij} } \mathbf{F}_2 +
                      \mathbf{F}_3 \right)                                \\[3ex]
              \mathbf{F}_k &=& && &&
                      \dfrac{d \sin \theta}{ | \mathbf{r}_\perp | } \mathbf{F}_2  \\[3ex]
              \end{array}                                             \\[5ex]
              \mbox{where ~}
              \mathbf{F}_1 = \mathbf{F}_{s} -
                        \dfrac{ \mathbf{r}_ij \cdot \mathbf{F}_{s} }
                              { \mathbf{r}_ij \cdot \mathbf{r}_{ij} } \mathbf{r}_{ij}
              \mbox{\,, ~}
              \mathbf{F}_2 = \mathbf{F}_1 -
                        \dfrac{ \mathbf{r}_\perp \cdot \mathbf{F}_{s} }
                              { \mathbf{r}_\perp \cdot \mathbf{r}_\perp } \mathbf{r}_\perp
              \mbox{~and ~}
              \mathbf{F}_3 = \dfrac{ \mathbf{r}_ij \cdot \mathbf{F}_{s} }
                               { \mathbf{r}_ij \cdot \mathbf{r}_{ij} } \mathbf{r}_\perp
      \end{array}

-  As a non-linear combination of three atoms, out of
   plane (:numref:`Fig. %s <fig-vsites>` 3out):

   .. math::

      \mathbf{r}_s ~=~ \mathbf{r}_i + a \mathbf{r}_ij + b \mathbf{r}_{ik} +
                      c (\mathbf{r}_ij \times \mathbf{r}_{ik})

-  This enables the construction of virtual sites out of the plane of
   the other atoms. The force on particles :math:`i,j` and :math:`k` due
   to the force on the virtual site can be computed as:

   .. math::

      \begin{array}{lcl}
      \mathbf{F}_j &=& \left[\begin{array}{ccc}
       a              &  -c\,z_{ik}   & c\,y_{ik}     \\[0.5ex]
       c\,z_{ik}      &   a           & -c\,x_{ik}    \\[0.5ex]
      -c\,y_{ik}      &   c\,x_{ik}   & a
      \end{array}\right]\mathbf{F}_{s}                                 \\
      \mathbf{F}_k &=& \left[\begin{array}{ccc}
       b              &   c\,z_{ij}   & -c\,y_{ij}    \\[0.5ex]
      -c\,z_{ij}      &   b           & c\,x_{ij}     \\[0.5ex]
       c\,y_{ij}      &  -c\,x_{ij}   & b
      \end{array}\right]\mathbf{F}_{s}                                 \\
      \mathbf{F}_i &=& \mathbf{F}_{s} - \mathbf{F}_j - \mathbf{F}_k
      \end{array}

-  From four atoms, with a fixed distance, see
   separate :numref:`Fig. %s <fig-vsite4fdn>`. This construction is a bit complex,
   in particular since the previous type (4fd) could be unstable which
   forced us to introduce a more elaborate construction:

.. _fig-vsite4fdn:

.. figure:: plots/vsite-4fdn.*
      :width: 5.00000cm

      The new 4fdn virtual site construction, which is stable even when
      all constructing atoms are in the same plane.

-
      .. math::   \begin{aligned}
                  \mathbf{r}_{ja} &=& a\, \mathbf{r}_{ik} - \mathbf{r}_{ij} = a\, (\mathbf{x}_k - \mathbf{x}_i) - (\mathbf{x}_j - \mathbf{x}_i) \nonumber \\
                  \mathbf{r}_{jb} &=& b\, \mathbf{r}_{il} - \mathbf{r}_{ij} = b\, (\mathbf{x}_l - \mathbf{x}_i) - (\mathbf{x}_j - \mathbf{x}_i) \nonumber \\
                  \mathbf{r}_m &=& \mathbf{r}_{ja} \times \mathbf{r}_{jb} \nonumber \\
                  \mathbf{x}_s &=& \mathbf{x}_i + c \frac{\mathbf{r}_m}{ | \mathbf{r}_m | }
                  \end{aligned}
                  :label: eqnvsite

-  In this case the virtual site is at a distance of :math:`|c|` from
   :math:`i`, while :math:`a` and :math:`b` are parameters. **Note**
   that the vectors :math:`\mathbf{r}_{ik}` and :math:`\mathbf{r}_{ij}`
   are not normalized to save floating-point operations. The force on
   particles :math:`i`, :math:`j`, :math:`k` and :math:`l` due to the
   force on the virtual site are computed through chain rule derivatives
   of the construction expression. This is exact and conserves energy,
   but it does lead to relatively lengthy expressions that we do not
   include here (over 200 floating-point operations). The interested
   reader can look at the source code in ``vsite.c``. Fortunately, this
   vsite type is normally only used for chiral centers such as
   :math:`C_{\alpha}` atoms in proteins.

   The new 4fdn construct is identified with a ‘type’ value of 2 in the
   topology. The earlier 4fd type is still supported internally (‘type’
   value 1), but it should not be used for new simulations. All current
   |Gromacs| tools will automatically generate type 4fdn instead.

-  A linear combination of :math:`N` atoms with relative
   weights :math:`a_i`. The weight for atom :math:`i` is:

   .. math:: w_i = a_i \left(\sum_{j=1}^N a_j \right)^{-1}

-   There are three options for setting the weights:

   -  center of geometry: equal weights

   -  center of mass: :math:`a_i` is the mass of atom :math:`i`; when in
      free-energy simulations the mass of the atom is changed, only the
      mass of the A-state is used for the weight

   -  center of weights: :math:`a_i` is defined by the user


.. _lrelstat:   

Long Range Electrostatics
-------------------------

Ewald summation
~~~~~~~~~~~~~~~

The total electrostatic energy of :math:`N` particles and their periodic
images is given by

.. math:: V=\frac{f}{2}\sum_{n_x}\sum_{n_y}
          \sum_{n_{z}*} \sum_{i}^{N} \sum_{j}^{N}
          \frac{q_i q_j}{{\bf r}_{ij,{\bf n}}}.
          :label: eqntotalcoulomb

:math:`(n_x,n_y,n_z)={\bf n}` is the box index vector, and the star
indicates that terms with :math:`i=j` should be omitted when
:math:`(n_x,n_y,n_z)=(0,0,0)`. The distance :math:`{\bf r}_{ij,{\bf n}}`
is the real distance between the charges and not the minimum-image. This
sum is conditionally convergent, but very slow.

Ewald summation was first introduced as a method to calculate long-range
interactions of the periodic images in crystals \ :ref:`105 <refEwald21>`. The idea
is to convert the single slowly-converging sum :eq:`eqn. %s <eqntotalcoulomb>`
into two quickly-converging terms and a constant term:

.. math::

   \begin{aligned}
   V &=& V_{\mathrm{dir}} + V_{\mathrm{rec}} + V_{0} \\[0.5ex]
   V_{\mathrm{dir}} &=& \frac{f}{2} \sum_{i,j}^{N}
   \sum_{n_x}\sum_{n_y}
   \sum_{n_{z}*} q_i q_j \frac{\mbox{erfc}(\beta {r}_{ij,{\bf n}} )}{{r}_{ij,{\bf n}}} \\[0.5ex]
   V_{\mathrm{rec}} &=& \frac{f}{2 \pi V} \sum_{i,j}^{N} q_i q_j
   \sum_{m_x}\sum_{m_y}
   \sum_{m_{z}*} \frac{\exp{\left( -(\pi {\bf m}/\beta)^2 + 2 \pi i
         {\bf m} \cdot ({\bf r}_i - {\bf r}_j)\right)}}{{\bf m}^2} \\[0.5ex]
   V_{0} &=& -\frac{f \beta}{\sqrt{\pi}}\sum_{i}^{N} q_i^2,\end{aligned}

where :math:`\beta` is a parameter that determines the relative weight
of the direct and reciprocal sums and :math:`{\bf m}=(m_x,m_y,m_z)`. In
this way we can use a short cut-off (of the order of :math:`1`  nm) in
the direct space sum and a short cut-off in the reciprocal space sum
(*e.g.* 10 wave vectors in each direction). Unfortunately, the
computational cost of the reciprocal part of the sum increases as
:math:`N^2` (or :math:`N^{3/2}` with a slightly better algorithm) and it
is therefore not realistic for use in large systems.

Using Ewald
^^^^^^^^^^^

Don’t use Ewald unless you are absolutely sure this is what you want -
for almost all cases the PME method below will perform much better. If
you still want to employ classical Ewald summation enter this in your
:ref:`mdp` file, if the side of your box is about :math:`3`  nm:

::

    coulombtype     = Ewald
    rvdw            = 0.9
    rlist           = 0.9
    rcoulomb        = 0.9
    fourierspacing  = 0.6
    ewald-rtol      = 1e-5

The ratio of the box dimensions and the fourierspacing parameter
determines the highest magnitude of wave vectors :math:`m_x,m_y,m_z` to
use in each direction. With a 3-nm cubic box this example would use
:math:`11` wave vectors (from :math:`-5` to :math:`5`) in each
direction. The ``ewald-rtol`` parameter is the relative strength of the
electrostatic interaction at the cut-off. Decreasing this gives you a
more accurate direct sum, but a less accurate reciprocal sum.

.. _pme:

PME
~~~

Particle-mesh Ewald is a method proposed by Tom
Darden \ :ref:`14 <refDarden93>` to improve the performance of the reciprocal sum.
Instead of directly summing wave vectors, the charges are assigned to a
grid using interpolation. The implementation in |Gromacs| uses cardinal
B-spline interpolation \ :ref:`15 <refEssmann95>`, which is referred to as
smooth PME (SPME). The grid is then Fourier transformed with a 3D FFT
algorithm and the reciprocal energy term obtained by a single sum over
the grid in k-space.

The potential at the grid points is calculated by inverse
transformation, and by using the interpolation factors we get the forces
on each atom.

The PME algorithm scales as :math:`N \log(N)`, and is substantially
faster than ordinary Ewald summation on medium to large systems. On very
small systems it might still be better to use Ewald to avoid the
overhead in setting up grids and transforms. For the parallelization of
PME see the section on MPMD PME (:ref:`mpmdpme`).

With the Verlet cut-off scheme, the PME direct space potential is
shifted by a constant such that the potential is zero at the cut-off.
This shift is small and since the net system charge is close to zero,
the total shift is very small, unlike in the case of the Lennard-Jones
potential where all shifts add up. We apply the shift anyhow, such that
the potential is the exact integral of the force.

Using PME
^^^^^^^^^

As an example for using Particle-mesh Ewald summation in |Gromacs|,
specify the following lines in your :ref:`mdp` file:

::

    coulombtype     = PME
    rvdw            = 0.9
    rlist           = 0.9
    rcoulomb        = 0.9
    fourierspacing  = 0.12
    pme-order       = 4
    ewald-rtol      = 1e-5

In this case the ``fourierspacing`` parameter determines the
maximum spacing for the FFT grid (i.e. minimum number of grid points),
and ``pme-order`` controls the interpolation order. Using
fourth-order (cubic) interpolation and this spacing should give
electrostatic energies accurate to about :math:`5\cdot10^{-3}`. Since
the Lennard-Jones energies are not this accurate it might even be
possible to increase this spacing slightly.

Pressure scaling works with PME, but be aware of the fact that
anisotropic scaling can introduce artificial ordering in some systems.

P3M-AD
~~~~~~

The Particle-Particle
Particle-Mesh
methods of Hockney & Eastwood can also be applied in |Gromacs| for the
treatment of long range electrostatic
interactions \ :ref:`106 <refHockney81>`. Although the P3M method was the first efficient long-range
electrostatics method for molecular simulation, the smooth PME (SPME)
method has largely replaced P3M as the method of choice in atomistic
simulations. One performance disadvantage of the original P3M method was
that it required 3 3D-FFT back transforms to obtain the forces on the
particles. But this is not required for P3M and the forces can be
derived through analytical differentiation of the potential, as done in
PME. The resulting method is termed P3M-AD. The only remaining
difference between P3M-AD and PME is the optimization of the lattice
Green influence function for error minimization that P3M uses. However,
in 2012 it has been shown that the SPME influence function can be
modified to obtain P3M \ :ref:`107 <refBallenegger2012>`. This means
that the advantage of error minimization in P3M-AD can be used at the
same computational cost and with the same code as PME, just by adding a
few lines to modify the influence function. However, at optimal
parameter setting the effect of error minimization in P3M-AD is less
than 10%. P3M-AD does show large accuracy gains with interlaced (also
known as staggered) grids, but that is not supported in |Gromacs| (yet).

P3M is used in |Gromacs| with exactly the same options as used with PME by
selecting the electrostatics type:

::

    coulombtype     = P3M-AD

Optimizing Fourier transforms and PME calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to optimize the parameters for calculation of
electrostatic interaction such as PME grid dimensions and cut-off radii.
This is particularly relevant to do before launching long production
runs.

:ref:`gmx mdrun` will automatically do a lot of PME
optimization, and |Gromacs| also includes a special tool,
:ref:`gmx tune_pme`, which automates the process of selecting
the optimal number of PME-only ranks.

Long Range Van der Waals interactions
-------------------------------------

Dispersion correction
~~~~~~~~~~~~~~~~~~~~~

In this section, we derive long-range corrections due to the use of a
cut-off for Lennard-Jones or Buckingham interactions. We assume that the
cut-off is so long that the repulsion term can safely be neglected, and
therefore only the dispersion term is taken into account. Due to the
nature of the dispersion interaction (we are truncating a potential
proportional to :math:`-r^{-6}`), energy and pressure corrections are
both negative. While the energy correction is usually small, it may be
important for free energy calculations where differences between two
different Hamiltonians are considered. In contrast, the pressure
correction is very large and can not be neglected under any
circumstances where a correct pressure is required, especially for any
NPT simulations. Although it is, in principle, possible to parameterize
a force field such that the pressure is close to the desired
experimental value without correction, such a method makes the
parameterization dependent on the cut-off and is therefore undesirable.

.. _ecorr:

Energy
^^^^^^

The long-range contribution of the dispersion interaction to the virial
can be derived analytically, if we assume a homogeneous system beyond
the cut-off distance :math:`r_c`. The dispersion energy between two
particles is written as:

.. math:: V({r_{ij}}) ~=~- C_6\,{r_{ij}}^{-6}

and the corresponding force is:

.. math:: \mathbf{F}_ij ~=~- 6\,C_6\,r_{ij}^{-8}\mathbf{r}_ij

In a periodic system it is not easy to calculate the full potentials,
so usually a cut-off is applied, which can be abrupt or smooth. We will
call the potential and force with cut-off :math:`V_c` and
:math:`\mathbf{F}_c`. The long-range contribution to the
dispersion energy in a system with :math:`N` particles and particle
density :math:`\rho` = :math:`N/V` is:

.. math:: V_{lr} ~=~ {\frac{1}{2}}N \rho\int_0^{\infty}   4\pi r^2 g(r) \left( V(r) -V_c(r) \right) {{{\rm d}r}}
          :label: eqnenercorr

We will integrate this for the shift function, which is the most
general form of van der Waals interaction available in |Gromacs|. The
shift function has a constant difference :math:`S` from 0 to :math:`r_1`
and is 0 beyond the cut-off distance :math:`r_c`. We can integrate
:eq:`eqn. %s <eqnenercorr>`, assuming that the density in the sphere within
:math:`r_1` is equal to the global density and the radial distribution
function :math:`g(r)` is 1 beyond :math:`r_1`:

.. math::

   \begin{aligned}
   \nonumber
   V_{lr}  &=& {\frac{1}{2}}N \left(
     \rho\int_0^{r_1}  4\pi r^2 g(r) \, C_6 \,S\,{{{\rm d}r}}
   + \rho\int_{r_1}^{r_c}  4\pi r^2 \left( V(r) -V_c(r) \right) {{{\rm d}r}}
   + \rho\int_{r_c}^{\infty}  4\pi r^2 V(r) \, {{{\rm d}r}}
   \right) \\
   & = & {\frac{1}{2}}N \left(\left(\frac{4}{3}\pi \rho r_1^{3} - 1\right) C_6 \,S
   + \rho\int_{r_1}^{r_c} 4\pi r^2 \left( V(r) -V_c(r) \right) {{{\rm d}r}}
   -\frac{4}{3} \pi N \rho\, C_6\,r_c^{-3}
   \right)\end{aligned}

where the term :math:`-1` corrects for the self-interaction. For a
plain cut-off we only need to assume that :math:`g(r)` is 1 beyond
:math:`r_c` and the correction reduces to \ :ref:`108 <refAllen87>`:

.. math::

   \begin{aligned}
   V_{lr} & = & -\frac{2}{3} \pi N \rho\, C_6\,r_c^{-3}\end{aligned}

If we consider, for example, a box of pure water, simulated with a
cut-off of 0.9 nm and a density of 1 g cm\ :math:`^{-3}` this correction
is :math:`-0.75` kJ mol\ :math:`^{-1}` per molecule.

For a homogeneous mixture we need to define an *average dispersion constant*:

.. math:: {\left< C_6 \right>}= \frac{2}{N(N-1)}\sum_i^N\sum_{j>i}^N C_6(i,j)\\
          :label: eqnavcsix

In |Gromacs|, excluded pairs of atoms do not contribute to the average.

In the case of inhomogeneous simulation systems, *e.g.* a system with a
lipid interface, the energy correction can be applied if
:math:`{\left< C_6 \right>}` for both components is comparable.

.. _virial:

Virial and pressure
^^^^^^^^^^^^^^^^^^^

The scalar virial of the system due to the dispersion interaction
between two particles :math:`i` and :math:`j` is given by:

.. math:: \Xi~=~-{\frac{1}{2}} \mathbf{r}_ij \cdot \mathbf{F}_ij ~=~ 3\,C_6\,r_{ij}^{-6}

The pressure is given by:

.. math:: P~=~\frac{2}{3\,V}\left(E_{kin} - \Xi\right)

The long-range correction to the virial is given by:

.. math:: \Xi_{lr} ~=~ {\frac{1}{2}}N \rho \int_0^{\infty} 4\pi r^2 g(r) (\Xi -\Xi_c) \,{{\rm d}r}

We can again integrate the long-range contribution to the virial
assuming :math:`g(r)` is 1 beyond :math:`r_1`:

.. math::

   \begin{aligned}
   \Xi_{lr}&=&	{\frac{1}{2}}N \rho \left(
       \int_{r_1}^{r_c}  4 \pi r^2 (\Xi -\Xi_c)  \,{{\rm d}r}+ \int_{r_c}^{\infty} 4 \pi r^2 3\,C_6\,{r_{ij}}^{-6}\,  {{\rm d}r}\right)	\nonumber\\
           &=&     {\frac{1}{2}}N \rho \left(
       \int_{r_1}^{r_c} 4 \pi r^2 (\Xi -\Xi_c) \, {{\rm d}r}+ 4 \pi C_6 \, r_c^{-3} \right)\end{aligned}

For a plain cut-off the correction to the pressure
is \ :ref:`108 <refAllen87>`:

.. math:: P_{lr}~=~-\frac{4}{3} \pi C_6\, \rho^2 r_c^{-3}

Using the same example of a water box, the correction to the virial is
0.75 kJ mol\ :math:`^{-1}` per molecule, the corresponding correction to
the pressure for SPC water is approximately :math:`-280` bar.

For homogeneous mixtures, we can again use the average dispersion
constant :math:`{\left< C_6 \right>}` (:eq:`eqn. %s <eqnavcsix>`):

.. math:: P_{lr}~=~-\frac{4}{3} \pi {\left< C_6 \right>}\rho^2 r_c^{-3}
          :label: eqnpcorr

For inhomogeneous systems, :eq:`eqn. %s <eqnpcorr>` can be applied under the
same restriction as holds for the energy (see sec. :ref:`ecorr`).

Lennard-Jones PME
~~~~~~~~~~~~~~~~~

In order to treat systems, using Lennard-Jones potentials, that are
non-homogeneous outside of the cut-off distance, we can instead use the
Particle-mesh Ewald method as discussed for electrostatics above. In
this case the modified Ewald equations become

.. math:: \begin{aligned}
          V &=& V_{\mathrm{dir}} + V_{\mathrm{rec}} + V_{0} \\[0.5ex]
          V_{\mathrm{dir}} &=& -\frac{1}{2} \sum_{i,j}^{N}
          \sum_{n_x}\sum_{n_y}
          \sum_{n_{z}*} \frac{C^{ij}_6 g(\beta {r}_{ij,{\bf n}})}{{r_{ij,{\bf n}}}^6}
          \end{aligned}
          :label: eqnljpmerealspace

.. math::
   \begin{aligned} 
   V_{\mathrm{rec}} &=& \frac{{\pi}^{\frac{3}{2}} \beta^{3}}{2V} \sum_{m_x}\sum_{m_y}\sum_{m_{z}*}
   f(\pi |{\mathbf m}|/\beta) \times \sum_{i,j}^{N} C^{ij}_6 {\mathrm{exp}}\left[-2\pi i {\bf m}\cdot({\bf r_i}-{\bf r_j})\right] \\[0.5ex]
   V_{0} &=& -\frac{\beta^{6}}{12}\sum_{i}^{N} C^{ii}_6\end{aligned}

where :math:`{\bf m}=(m_x,m_y,m_z)`, :math:`\beta` is the parameter
determining the weight between direct and reciprocal space, and
:math:`{C^{ij}_6}` is the combined dispersion parameter for particle
:math:`i` and :math:`j`. The star indicates that terms with
:math:`i = j` should be omitted when :math:`((n_x,n_y,n_z)=(0,0,0))`,
and :math:`{\bf r}_{ij,{\bf n}}` is the real distance between the
particles. Following the derivation by
Essmann \ :ref:`15 <refEssmann95>`, the functions :math:`f` and :math:`g`
introduced above are defined as

.. math::

   \begin{aligned}
   f(x)&=&1/3\left[(1-2x^2){\mathrm{exp}}(-x^2) + 2{x^3}\sqrt{\pi}\,{\mathrm{erfc}}(x) \right] \\
   g(x)&=&{\mathrm{exp}}(-x^2)(1+x^2+\frac{x^4}{2}).\end{aligned}

The above methodology works fine as long as the dispersion parameters
can be combined geometrically (:eq:`eqn. %s <eqncomb>`) in the same way as the
charges for electrostatics

.. math:: C^{ij}_{6,\mathrm{geom}} = \left(C^{ii}_6 \, C^{jj}_6\right)^{1/2}

For Lorentz-Berthelot combination rules (:eq:`eqn. %s <eqnlorentzberthelot>`),
the reciprocal part of this sum has to be calculated seven times due to
the splitting of the dispersion parameter according to

.. math:: C^{ij}_{6,\mathrm{L-B}} = (\sigma_i+\sigma_j)^6=\sum_{n=0}^{6} P_{n}\sigma_{i}^{n}\sigma_{j}^{(6-n)},

for :math:`P_{n}` the Pascal triangle coefficients. This introduces a
non-negligible cost to the reciprocal part, requiring seven separate
FFTs, and therefore this has been the limiting factor in previous
attempts to implement LJ-PME. A solution to this problem is to use
geometrical combination rules in order to calculate an approximate
interaction parameter for the reciprocal part of the potential, yielding
a total interaction of

.. math::

   \begin{aligned}
   V(r<r_c) & = & \underbrace{C^{\mathrm{dir}}_6 g(\beta r) r^{-6}}_{\mathrm{Direct \  space}} + \underbrace{C^\mathrm{recip}_{6,\mathrm{geom}} [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}} \nonumber \\
   &=& C^\mathrm{recip}_{6,\mathrm{geom}}r^{-6} + \left(C^{\mathrm{dir}}_6-C^\mathrm{recip}_{6,\mathrm{geom}}\right)g(\beta r)r^{-6} \\
   V(r>r_c) & = & \underbrace{C^\mathrm{recip}_{6,\mathrm{geom}} [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}}.\end{aligned}

This will preserve a well-defined Hamiltonian and significantly
increase the performance of the simulations. The approximation does
introduce some errors, but since the difference is located in the
interactions calculated in reciprocal space, the effect will be very
small compared to the total interaction energy. In a simulation of a
lipid bilayer, using a cut-off of 1.0 nm, the relative error in total
dispersion energy was below 0.5%. A more thorough discussion of this can
be found in :ref:`109 <refWennberg13>`.

In |Gromacs| we now perform the proper calculation of this interaction by
subtracting, from the direct-space interactions, the contribution made
by the approximate potential that is used in the reciprocal part

.. math:: V_\mathrm{dir} = C^{\mathrm{dir}}_6 r^{-6} - C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}.
          :label: eqnljpmedirectspace

This potential will reduce to the expression in
:eq:`eqn. %s <eqnljpmerealspace>` when
:math:`C^{\mathrm{dir}}_6 = C^\mathrm{recip}_6`, and the total
interaction is given by

.. math:: \begin{aligned}
          \nonumber V(r<r_c) &=& \underbrace{C^{\mathrm{dir}}_6 r^{-6} - C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}}_{\mathrm{Direct \  space}} + \underbrace{C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}} \\ 
          &=&C^{\mathrm{dir}}_6 r^{-6}
          \end{aligned}
          :label: eqnljpmecorr2

.. math::

   \begin{aligned} 
   V(r>r_c) &=& C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}.\end{aligned}

For the case when :math:`C^{\mathrm{dir}}_6 \neq C^\mathrm{recip}_6`
this will retain an unmodified LJ force up to the cut-off, and the error
is an order of magnitude smaller than in simulations where the
direct-space interactions do not account for the approximation used in
reciprocal space. When using a VdW interaction modifier of
potential-shift, the constant

.. math:: \left(-C^{\mathrm{dir}}_6 + C^\mathrm{recip}_6 [1 - g(\beta r_c)]\right) r_c^{-6}

is added to :eq:`eqn. %s <eqnljpmecorr2>` in order to ensure that the potential
is continuous at the cutoff. Note that, in the same way as
:eq:`eqn. %s <eqnljpmedirectspace>`, this degenerates into the expected
:math:`-C_6g(\beta r_c)r^{-6}_c` when :math:`C^{\mathrm{dir}}_6 =
C^\mathrm{recip}_6`. In addition to this, a long-range dispersion
correction can be applied to correct for the approximation using a
combination rule in reciprocal space. This correction assumes, as for
the cut-off LJ potential, a uniform particle distribution. But since the
error of the combination rule approximation is very small this
long-range correction is not necessary in most cases. Also note that
this homogenous correction does not correct the surface tension, which
is an inhomogeneous property.

Using LJ-PME
^^^^^^^^^^^^

As an example for using Particle-mesh Ewald summation for Lennard-Jones
interactions in |Gromacs|, specify the following lines in your :ref:`mdp` file:

::

    vdwtype          = PME
    rvdw             = 0.9
    vdw-modifier     = Potential-Shift
    rlist            = 0.9
    rcoulomb         = 0.9
    fourierspacing   = 0.12
    pme-order        = 4
    ewald-rtol-lj    = 0.001
    lj-pme-comb-rule = geometric

The same Fourier grid and interpolation order are used if both LJ-PME
and electrostatic PME are active, so the settings for
``fourierspacing`` and ``pme-order`` are common
to both. ``ewald-rtol-lj`` controls the splitting between
direct and reciprocal space in the same way as
``ewald-rtol``. In addition to this, the combination rule to
be used in reciprocal space is determined by
``lj-pme-comb-rule``. If the current force field uses
Lorentz-Berthelot combination rules, it is possible to set
``lj-pme-comb-rule = geometric`` in order to gain a
significant increase in performance for a small loss in accuracy. The
details of this approximation can be found in the section above.

Note that the use of a complete long-range dispersion correction means
that as with Coulomb PME, ``rvdw`` is now a free parameter
in the method, rather than being necessarily restricted by the
force-field parameterization scheme. Thus it is now possible to optimize
the cutoff, spacing, order and tolerance terms for accuracy and best
performance.

Naturally, the use of LJ-PME rather than LJ cut-off adds computation and
communication done for the reciprocal-space part, so for best
performance in balancing the load of parallel simulations using PME-only
ranks, more such ranks should be used. It may be possible to improve
upon the automatic load-balancing used by :ref:`mdrun <gmx mdrun>`.

Force field
-----------

A force field is built up from two distinct components:

-  The set of equations (called the *potential functions*) used to
   generate the potential energies and their derivatives, the forces.
   These are described in detail in the previous chapter.

-  The parameters used in this set of equations. These are not given in
   this manual, but in the data files corresponding to your |Gromacs|
   distribution.

Within one set of equations various sets of parameters can be used. Care
must be taken that the combination of equations and parameters form a
consistent set. It is in general dangerous to make *ad hoc* changes in a
subset of parameters, because the various contributions to the total
force are usually interdependent. This means in principle that every
change should be documented, verified by comparison to experimental data
and published in a peer-reviewed journal before it can be used.

|Gromacs| |version| includes several force fields, and
additional ones are available on the website. If you do not know which
one to select we recommend GROMOS-96 for united-atom setups and
OPLS-AA/L for all-atom parameters. That said, we describe the available
options in some detail.

All-hydrogen force field
~~~~~~~~~~~~~~~~~~~~~~~~

The GROMOS-87-based all-hydrogen force field is almost identical to the
normal GROMOS-87 force field, since the extra hydrogens have no
Lennard-Jones interaction and zero charge. The only differences are in
the bond angle and improper dihedral angle terms. This force field is
only useful when you need the exact hydrogen positions, for instance for
distance restraints derived from NMR measurements. When citing this
force field please read the previous paragraph.

GROMOS-96
~~~~~~~~~

|Gromacs| supports the GROMOS-96 force fields \ :ref:`77 <refgromos96>`. All
parameters for the 43A1, 43A2 (development, improved alkane dihedrals),
45A3, 53A5, and 53A6 parameter sets are included. All standard building
blocks are included and topologies can be built automatically by
:ref:`pdb2gmx <gmx pdb2gmx>`.

The GROMOS-96 force field is a further development of the GROMOS-87
force field. It has improvements over the GROMOS-87 force field for
proteins and small molecules. **Note** that the sugar parameters present
in 53A6 do correspond to those published in 2004\ :ref:`110 <refOostenbrink2004>`,
which are different from those present in 45A4, which is not
included in |Gromacs| at this time. The 45A4 parameter set corresponds to
a later revision of these parameters. The GROMOS-96 force field is not,
however, recommended for use with long alkanes and lipids. The GROMOS-96
force field differs from the GROMOS-87 force field in a few respects:

-  the force field parameters

-  the parameters for the bonded interactions are not linked to atom
   types

-  a fourth power bond stretching potential (:ref:`G96bond`)

-  an angle potential based on the cosine of the angle
   (:ref:`G96angle`)

There are two differences in implementation between |Gromacs| and
GROMOS-96 which can lead to slightly different results when simulating
the same system with both packages:

-  in GROMOS-96 neighbor searching for solvents is performed on the
   first atom of the solvent molecule. This is not implemented in
   |Gromacs|, but the difference with searching by centers of charge
   groups is very small

-  the virial in GROMOS-96 is molecule-based. This is not implemented in
   |Gromacs|, which uses atomic virials

The GROMOS-96 force field was parameterized with a Lennard-Jones cut-off
of 1.4 nm, so be sure to use a Lennard-Jones cut-off
(``rvdw``) of at least 1.4. A larger cut-off is possible
because the Lennard-Jones potential and forces are almost zero beyond
1.4 nm.

GROMOS-96 files
^^^^^^^^^^^^^^^

|Gromacs| can read and write GROMOS-96 coordinate and trajectory files.
These files should have the extension :ref:`g96`. Such a file
can be a GROMOS-96 initial/final configuration file, a coordinate
trajectory file, or a combination of both. The file is fixed format; all
floats are written as 15.9, and as such, files can get huge. |Gromacs|
supports the following data blocks in the given order:

-  Header block:

   ::

       TITLE (mandatory)

-  Frame blocks:

   ::

       TIMESTEP (optional)
       POSITION/POSITIONRED (mandatory)
       VELOCITY/VELOCITYRED (optional)
       BOX (optional)

See the GROMOS-96 manual \ :ref:`77 <refgromos96>` for a complete
description of the blocks. **Note** that all |Gromacs| programs can read
compressed (.Z) or gzipped (.gz) files.

OPLS/AA
~~~~~~~

AMBER
~~~~~

|Gromacs| provides native support for the following AMBER force fields:

-  AMBER94 \ :ref:`111 <refCornell1995>`

-  AMBER96 \ :ref:`112 <refKollman1996>`

-  AMBER99 \ :ref:`113 <refWang2000>`

-  AMBER99SB \ :ref:`114 <refHornak2006>`

-  AMBER99SB-ILDN \ :ref:`115 <refLindorff2010>`

-  AMBER03 \ :ref:`116 <refDuan2003>`

-  AMBERGS \ :ref:`117 <refGarcia2002>`

.. _charmmff:

CHARMM
~~~~~~

|Gromacs| supports the CHARMM force field for
proteins \ :ref:`118 <refmackerell04>`, :ref:`119 <refmackerell98>`,
lipids \ :ref:`120 <reffeller00>` and nucleic
acids \ :ref:`121 <reffoloppe00>`, :ref:`122 <refMac2000>`. The protein
parameters (and to some extent
the lipid and nucleic acid parameters) were thoroughly tested – both by
comparing potential energies between the port and the standard parameter
set in the CHARMM molecular simulation package, as well by how the
protein force field behaves together with |Gromacs|-specific techniques
such as virtual sites (enabling long time steps) recently
implemented \ :ref:`123 <refLarsson10>` – and the details and results are
presented in the paper by Bjelkmar et al. \ :ref:`124 <refBjelkmar10>`.
The nucleic acid parameters, as well as the ones for HEME, were
converted and tested by Michel Cuendet.

When selecting the CHARMM force field in
:ref:`pdb2gmx <gmx pdb2gmx>` the default option
is to use CMAP (for torsional correction map).
To exclude CMAP, use ``-nocmap``. The basic form of the CMAP
term implemented in |Gromacs| is a function of the :math:`\phi` and
:math:`\psi` backbone torsion angles. This term is defined in the
``rtp`` file by a ``[ cmap ]`` statement at the
end of each residue supporting CMAP. The following five atom names
define the two torsional angles. Atoms 1-4 define :math:`\phi`, and
atoms 2-5 define :math:`\psi`. The corresponding atom types are then
matched to the correct CMAP type in the ``cmap.itp`` file
that contains the correction maps.

A port of the CHARMM36 force field for use with |Gromacs| is also
available at `the MacKerell lab webpage <http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs>`_.

For branched polymers or other topologies not supported by
:ref:`pdb2gmx <gmx pdb2gmx>`, it is possible to
use TopoTools \ :ref:`125 <refkohlmeyer2016>` to generate a |Gromacs| top
file.

.. _cgforcefields:

Coarse-grained force fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coarse-graining is a systematic way of reducing the
number of degrees of freedom representing a system of interest. To
achieve this, typically whole groups of atoms are represented by single
beads and the coarse-grained force fields describes their effective
interactions. Depending on the choice of parameterization, the
functional form of such an interaction can be complicated and often
tabulated potentials are used.

Coarse-grained models are designed to reproduce certain properties of a
reference system. This can be either a full atomistic model or even
experimental data. Depending on the properties to reproduce there are
different methods to derive such force fields. An incomplete list of
methods is given below:

-  Conserving free energies

   -  Simplex method

   -  MARTINI force field (see next section)

-  Conserving distributions (like the radial distribution function),
   so-called structure-based coarse-graining

   -  (iterative) Boltzmann inversion

   -  Inverse Monte Carlo

-  Conversing forces

   -  Force matching

Note that coarse-grained potentials are state dependent (e.g.
temperature, density,...) and should be re-parametrized depending on the
system of interest and the simulation conditions. This can for example
be done using the Versatile Object-oriented Toolkit for Coarse-Graining
Applications (VOTCA) (**???**). The package was designed to assists in
systematic coarse-graining, provides implementations for most of the
algorithms mentioned above and has a well tested interface to |Gromacs|.
It is available as open source and further information can be found at
`www.votca.org <http://www.votca.org>`_.

MARTINI
~~~~~~~

The MARTINI force field is a coarse-grain parameter set that allows for
the construction of many systems, including proteins and membranes.

PLUM
~~~~

The PLUM force field :ref:`126 <refbereau12>` is an example of a solvent-free
protein-membrane model for which the membrane was derived from
structure-based coarse-graining \ :ref:`127 <refwang_jpcb10>`. A |Gromacs|
implementation can be found at
`code.google.com/p/plumx <http://code.google.com/p/plumx/>`__.

