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
          :label: eqnharmbondstretch

See also :numref:`Fig. %s <fig-bstretch1>`, with the force given by:

.. math:: \mathbf{F}_i(\mathbf{r}_{ij}) = k^b_{ij}({r_{ij}}-b_{ij}) {\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}
          :label: eqnharmbondstretchforce

.. _g96bond:

Fourth power potential
^^^^^^^^^^^^^^^^^^^^^^

In the GROMOS-96 force field \ :ref:`77 <refgromos96>`, the covalent bond
potential is, for reasons of computational efficiency, written as:

.. math:: V_b~({r_{ij}}) = \frac{1}{4}k^b_{ij}\left({r_{ij}}^2-b_{ij}^2\right)^2
          :label: eqng96bond

The corresponding force is:

.. math:: \mathbf{F}_i(\mathbf{r}_{ij}) = k^b_{ij}({r_{ij}}^2-b_{ij}^2)~\mathbf{r}_ij
          :label: eqng96bondforce

The force constants for this form of the potential are related to the
usual harmonic force constant :math:`k^{b,\mathrm{harm}}`
(sec. :ref:`bondpot`) as

.. math:: 2 k^b b_{ij}^2 = k^{b,\mathrm{harm}}
          :label: eqn96harmrelation

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
          :label: eqnmorsebond

See also :numref:`Fig. %s <fig-morse>`, and the corresponding force is:

.. math:: \begin{array}{rcl}
          \displaystyle {\bf F}_{morse} ({\bf r}_{ij})&=&2 D_{ij} \beta_{ij} \exp(-\beta_{ij}(r_{ij}-b_{ij})) * \\
          \displaystyle \: & \: &[1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))] \frac{\displaystyle {\bf r}_{ij}}{\displaystyle r_{ij}},
          \end{array}
          :label: eqnmorsebondforce

where :math:`\displaystyle D_{ij}`  is the depth of the well in
kJ/mol, :math:`\displaystyle \beta_{ij}` defines the steepness of the
well (in nm\ :math:`^{-1}`), and :math:`\displaystyle b_{ij}` is the
equilibrium distance in nm. The steepness parameter
:math:`\displaystyle \beta_{ij}` can be expressed in terms of the reduced mass of the atoms *i* and
*j*, the fundamental vibration frequency :math:`\displaystyle\omega_{ij}` and the well depth :math:`\displaystyle D_{ij}`:

.. math:: \displaystyle \beta_{ij}= \omega_{ij} \sqrt{\frac{\mu_{ij}}{2 D_{ij}}}
          :label: eqnmorsefreq

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

.. math:: \begin{array}{rcl}
          \displaystyle V_{morse} (r_{ij})&=&D_{ij} [1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))]^2\\
          \displaystyle \:&=&D_{ij} [1 - (1 -\sqrt{\frac{k_{ij}}{2 D_{ij}}}(r_{ij}-b_{ij}))]^2\\
          \displaystyle \:&=&\frac{1}{2} k_{ij} (r_{ij}-b_{ij}))^2
          \end{array}
          :label: eqnharmfrommorse

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
          :label: eqncubicbond

A flexible water model (based on the SPC water model \ :ref:`80 <refBerendsen81>`)
including a cubic bond stretching potential for the O-H bond was
developed by Ferguson \ :ref:`81 <refFerguson95>`. This model was found to yield a
reasonable infrared spectrum. The Ferguson water model is available in
the |Gromacs| library (``flexwat-ferguson.itp``). It should be noted that the
potential is asymmetric: overstretching leads to infinitely low
energies. The integration timestep is therefore limited to 1 fs.

The force corresponding to this potential is:

.. math:: \mathbf{F}_i(\mathbf{r}_{ij}) = 2k^b_{ij}({r_{ij}}-b_{ij})~{\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}+ 3k^b_{ij}k^{cub}_{ij}({r_{ij}}-b_{ij})^2~{\frac{{\mathbf{r}_{ij}}}{{r_{ij}}}}
          :label: eqncubicbondforce

FENE bond stretching potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In coarse-grained polymer simulations the beads are often connected by a
FENE (finitely extensible nonlinear elastic) potential \ :ref:`82 <refWarner72>`:

.. math:: V_{\mbox{FENE}}({r_{ij}}) =
          -{\frac{1}{2}}k^b_{ij} b^2_{ij} \log\left(1 - \frac{{r_{ij}}^2}{b^2_{ij}}\right)
          :label: eqnfenebond

The potential looks complicated, but the expression for the force is
simpler:

.. math:: F_{\mbox{FENE}}(\mathbf{r}_{ij}) =
          -k^b_{ij} \left(1 - \frac{{r_{ij}}^2}{b^2_{ij}}\right)^{-1} \mathbf{r}_{ij}
          :label: eqnfenebondforce

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

   Principle of angle vibration (left) and the bond angle potential.

.. math:: V_a({\theta_{ijk}}) = {\frac{1}{2}}k^{\theta}_{ijk}({\theta_{ijk}}-{\theta_{ijk}}^0)^2
          :label: eqnharmangle

As the bond-angle vibration is represented by a harmonic potential, the
form is the same as the bond stretching
(:numref:`Fig. %s <fig-bstretch1>`).

The force equations are given by the chain rule:

.. math:: \begin{array}{l}
          \mathbf{F}_i    ~=~ -\displaystyle\frac{d V_a({\theta_{ijk}})}{d \mathbf{r}_i}   \\
          \mathbf{F}_k    ~=~ -\displaystyle\frac{d V_a({\theta_{ijk}})}{d \mathbf{r}_k}   \\
          \mathbf{F}_j    ~=~ -\mathbf{F}_i-\mathbf{F}_k
          \end{array}
          ~ \mbox{ ~ where ~ } ~
           {\theta_{ijk}}= \arccos \frac{(\mathbf{r}_{ij} \cdot \mathbf{r}_{kj})}{r_{ij}r_{kj}}
          :label: eqnharmangleforce

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
          :label: eqnG96angle

where

.. math:: \cos({\theta_{ijk}}) = \frac{\mathbf{r}_{ij}\cdot\mathbf{r}_{kj}}{{r_{ij}}r_{kj}}
          :label: eqnG96anglecos

The corresponding force can be derived by partial differentiation with
respect to the atomic positions. The force constants in this function
are related to the force constants in the harmonic form
:math:`k^{\theta,\mathrm{harm}}` (:ref:`harmonicangle`) by:

.. math:: k^{\theta} \sin^2({\theta_{ijk}}^0) = k^{\theta,\mathrm{harm}}
          :label: eqnG96angleFC

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
:math:`180^{\circ}` value, the bending potential :eq:`eqn %s <eqnG96angle>` is
divided by a :math:`\sin^2\theta` factor:

.. math:: V_{\rm ReB}(\theta_i) = \frac{1}{2} k_{\theta} \frac{(\cos\theta_i - \cos\theta_0)^2}{\sin^2\theta_i}.
          :label: eqnReB

:numref:`Figure %s <fig-ReB>` shows the comparison between the ReB potential,
:eq:`%s <eqnReB>`, and the standard one :eq:`%s <eqnG96angle>`.

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

.. math:: F_{\rm ReB}(\theta_i) = -\frac{k_{\theta}}{\sin^4\theta_i}(\cos\theta_i - \cos\theta_0) (1 - \cos\theta_i\cos\theta_0) \frac{\partial \cos\theta_i}{\partial \vec r_{k}}.
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
          :label: eqnUBAngle

The force equations can be deduced from sections :ref:`harmonicbond`
and :ref:`harmonicangle`.

Linear Angle potential
~~~~~~~~~~~~~~~~~~~~~~

The linear angle potential was designed especially for linear compounds
such as nitriles and for carbon dioxide \ :ref:`190 <refSpoel2020>`. 
It avoids the calculation of the angle per se, since the angle force
is not well-defined if the angle is 180 degrees. Rather, it computes the
deviation of a central atom in a triplet *i,j,k* from a reference position

.. math:: \mathbf{x}_j^0 = a \mathbf{x}_i + (1-a) \mathbf{x}_k

where a is defined by the bond-length *i-j* and *j-k*, in a symmetric
molecule such as carbon dioxide *a = 1/2*. If the compound has different
bond lengths :math:`b_{ij}` and :math:`b_{jk}` respectivey, we instead have

.. math:: a = \frac{b_{jk}}{b_{ij}+b_{jk}}.

If the order of atoms is changed to *k,j,i*, *a* needs to be 
replaced by *1-a*. The energy is now given by

.. math:: V_{lin} = \frac{k_{lin}}{2}\left(\mathbf{x}_j - \mathbf{x}_j^0\right)^2

with :math:`k_{lin}` the force constant. For examples, and a derivation of the forces from the energy function, see ref. \ :ref:`190 <refSpoel2020>`. 

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
          :label: eqncrossbbforce

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
          :label: eqncrossba

where :math:`k_{r\theta}` is the force constant, :math:`r_{3e}` is the
:math:`i-k` distance, and the other constants are the same as in
:eq:`Equation %s <eqncrossbb>`. The force associated with the potential on atom
:math:`i` is:

.. math:: \mathbf{F}_{i} ~=~ -k_{r\theta}
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
          :label: eqncrossbaforce

Quartic angle potential
~~~~~~~~~~~~~~~~~~~~~~~

For special purposes there is an angle potential that uses a fourth
order polynomial:

.. math:: V_q({\theta_{ijk}}) ~=~ \sum_{n=0}^5 C_n ({\theta_{ijk}}-{\theta_{ijk}}^0)^n
          :label: eqnquarticangle

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
          :label: eqnharmimpdihedral

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
          :label: eqnperiodicpropdihedral

Proper dihedrals: Ryckaert-Bellemans function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| For alkanes, the following proper dihedral potential is often used
  (see :numref:`Fig. %s <fig-rbdih>`):

  .. math:: V_{rb}(\phi_{ijkl}) = \sum_{n=0}^5 C_n( \cos(\psi ))^n,
            :label: eqnRBproperdihedral

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

  .. math:: V_{rb} (\phi_{ijkl}) ~=~ \frac{1}{2} \left[F_1(1+\cos(\phi)) + F_2(
            1-\cos(2\phi)) + F_3(1+\cos(3\phi)) + F_4(1-\cos(4\phi))\right]
            :label: eqnRBproperdihedralFourier

| Because of the equalities :math:`\cos(2\phi) = 2\cos^2(\phi) - 1`,
  :math:`\cos(3\phi) = 4\cos^3(\phi) - 3\cos(\phi)` and
  :math:`\cos(4\phi) = 8\cos^4(\phi) - 8\cos^2(\phi) + 1` one can
  translate the OPLS parameters to Ryckaert-Bellemans parameters as
  follows:

  .. math:: \displaystyle
            \begin{array}{rcl}
            \displaystyle C_0&=&F_2 + \frac{1}{2} (F_1 + F_3)\\
            \displaystyle C_1&=&\frac{1}{2} (- F_1 + 3 \, F_3)\\
            \displaystyle C_2&=& -F_2 + 4 \, F_4\\
            \displaystyle C_3&=&-2 \, F_3\\
            \displaystyle C_4&=&-4 \, F_4\\
            \displaystyle C_5&=&0
            \end{array}
            :label: eqnoplsRBconversion

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

  .. math:: V_{F} (\phi_{ijkl}) ~=~ \frac{1}{2} \left[C_1(1+\cos(\phi)) + C_2(
            1-\cos(2\phi)) + C_3(1+\cos(3\phi)) + C_4(1-\cos(4\phi))\right],
            :label: eqnfourierproperdihedral

| Internally, |Gromacs| uses the Ryckaert-Bellemans code to compute
  Fourier dihedrals (see above), because this is more efficient.
| **Note:** Mind the conversion from **kcal mol**\ :math:`^{-1}` for
  literature OPLS parameters to **kJ mol**\ :math:`^{-1}` in |Gromacs|.

Proper dihedrals: Restricted torsion potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a manner very similar to the restricted bending potential (see
:ref:`ReB`), a restricted torsion/dihedral potential is introduced:

.. math:: V_{\rm ReT}(\phi_i) = \frac{1}{2} k_{\phi} \frac{(\cos\phi_i - \cos\phi_0)^2}{\sin^2\phi_i}
          :label: eqnReT

with the advantages of being a function of :math:`\cos\phi` (no
problems taking the derivative of :math:`\sin\phi`) and of keeping the
torsion angle at only one minimum value. In this case, the factor
:math:`\sin^2\phi` does not allow the dihedral angle to move from the
[:math:`-180^{\circ}`:0] to [0::math:`180^{\circ}`] interval, i.e. it
cannot have maxima both at :math:`-\phi_0` and :math:`+\phi_0` maxima,
but only one of them. For this reason, all the dihedral angles of the
starting configuration should have their values in the desired angles
interval and the equilibrium :math:`\phi_0` value should not be too
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
          :label: eqnCBT

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
   (:eq:`%s <eqnCBT>` with :math:`k = 10` kJ mol\ :math:`^{-1}`,
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
          :label: eqnforcecbt

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

Bonded pair and 1-4 interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most force fields do not use normal Lennard-Jones and Coulomb interactions
for atoms separated by three bonds, the so-called 1-4 interactions. These
interactions are still affected by the modified electronic distributions
due to the chemical bonds and they are modified in the force field by
the dihedral terms. For this reason the Lennard-Jones and Coulomb 1-4
interactions are often scaled down, by a fixed factor given by the force field.
These factors can be supplied in the topology and the parameters can also
be overriden per 1-4 interaction or atom type pair. The pair interactions
can be used for any atom pair in a molecule, not only 1-4 pairs.
The non-bonded interactions between such pairs should be excluded to avoid double
interactions. Plain Lennard-Jones and Coulomb interactions are used which
are not affected by the non-bonded interaction treatment and potential
modifiers.

Tabulated bonded interaction functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| For full flexibility, any functional shape can be used for bonds,
  angles and dihedrals through user-supplied tabulated functions. The
  functional shapes are:

  .. math:: \begin{aligned}
            V_b(r_{ij})      &=& k \, f^b_n(r_{ij}) \\
            V_a({\theta_{ijk}})       &=& k \, f^a_n({\theta_{ijk}}) \\
            V_d(\phi_{ijkl}) &=& k \, f^d_n(\phi_{ijkl})\end{aligned}
            :label: eqntabuöatedbond

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
