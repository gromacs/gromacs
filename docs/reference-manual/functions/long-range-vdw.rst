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
          :label: eqnlrljEdisp

and the corresponding force is:

.. math:: \mathbf{F}_{ij} ~=~- 6\,C_6\,r_{ij}^{-8}\mathbf{r}_{ij}
          :label: eqnlrljFdisp

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

.. math:: \begin{aligned}
          V_{lr}  &=& {\frac{1}{2}}N \left(
            \rho\int_0^{r_1}  4\pi r^2 g(r) \, C_6 \,S\,{{{\rm d}r}}
          + \rho\int_{r_1}^{r_c}  4\pi r^2 \left( V(r) -V_c(r) \right) {{{\rm d}r}}
          + \rho\int_{r_c}^{\infty}  4\pi r^2 V(r) \, {{{\rm d}r}}
          \right) \\
          & = & {\frac{1}{2}}N \left(\left(\frac{4}{3}\pi \rho r_1^{3} - 1\right) C_6 \,S
          + \rho\int_{r_1}^{r_c} 4\pi r^2 \left( V(r) -V_c(r) \right) {{{\rm d}r}}
          -\frac{4}{3} \pi N \rho\, C_6\,r_c^{-3}
          \right)\end{aligned}
          :label: eqnlrljshift

where the term :math:`-1` corrects for the self-interaction. For a
plain cut-off we only need to assume that :math:`g(r)` is 1 beyond
:math:`r_c` and the correction reduces to \ :ref:`108 <refAllen87>`:

.. math:: \begin{aligned}
          V_{lr} & = & -\frac{2}{3} \pi N \rho\, C_6\,r_c^{-3}\end{aligned}
          :label: eqnlrljcorrreduced

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

.. math:: \Xi~=~-{\frac{1}{2}} \mathbf{r}_{ij} \cdot \mathbf{F}_{ij} ~=~ 3\,C_6\,r_{ij}^{-6}
          :label: eqnlrljdispvirial

The pressure is given by:

.. math:: P~=~\frac{2}{3\,V}\left(E_{kin} - \Xi\right)
          :label: eqnlrljpressure

The long-range correction to the virial is given by:

.. math:: \Xi_{lr} ~=~ {\frac{1}{2}}N \rho \int_0^{\infty} 4\pi r^2 g(r) (\Xi -\Xi_c) \,{{\rm d}r}
          :label: eqnlrljcorrvirial

We can again integrate the long-range contribution to the virial
assuming :math:`g(r)` is 1 beyond :math:`r_1`:

.. math:: \begin{aligned}
          \Xi_{lr}&=&	{\frac{1}{2}}N \rho \left(
              \int_{r_1}^{r_c}  4 \pi r^2 (\Xi -\Xi_c)  \,{{\rm d}r}+ \int_{r_c}^{\infty} 4 \pi r^2 3\,C_6\,{r_{ij}}^{-6}\,  {{\rm d}r}\right)	\nonumber\\
                  &=&     {\frac{1}{2}}N \rho \left(
              \int_{r_1}^{r_c} 4 \pi r^2 (\Xi -\Xi_c) \, {{\rm d}r}+ 4 \pi C_6 \, r_c^{-3} \right)\end{aligned}
          :label: eqnlrljvirialcontrib

For a plain cut-off the correction to the pressure
is \ :ref:`108 <refAllen87>`:

.. math:: P_{lr}~=~-\frac{4}{3} \pi C_6\, \rho^2 r_c^{-3}
          :label: eqnlrljpressurecorr

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

.. math:: \begin{aligned} 
          V_{\mathrm{rec}} &=& \frac{{\pi}^{\frac{3}{2}} \beta^{3}}{2V} \sum_{m_x}\sum_{m_y}\sum_{m_{z}*}
          f(\pi | {\mathbf m} | /\beta) \times \sum_{i,j}^{N} C^{ij}_6 {\mathrm{exp}}\left[-2\pi i {\bf m}\cdot({\bf r_i}-{\bf r_j})\right] \\[0.5ex]
          V_{0} &=& -\frac{\beta^{6}}{12}\sum_{i}^{N} C^{ii}_6\end{aligned}
          :label: eqnljpmerealspace2

where :math:`{\bf m}=(m_x,m_y,m_z)`, :math:`\beta` is the parameter
determining the weight between direct and reciprocal space, and
:math:`{C^{ij}_6}` is the combined dispersion parameter for particle
:math:`i` and :math:`j`. The star indicates that terms with
:math:`i = j` should be omitted when :math:`((n_x,n_y,n_z)=(0,0,0))`,
and :math:`{\bf r}_{ij,{\bf n}}` is the real distance between the
particles. Following the derivation by
Essmann \ :ref:`15 <refEssmann95>`, the functions :math:`f` and :math:`g`
introduced above are defined as

.. math:: \begin{aligned}
          f(x)&=&1/3\left[(1-2x^2){\mathrm{exp}}(-x^2) + 2{x^3}\sqrt{\pi}\,{\mathrm{erfc}}(x) \right] \\
          g(x)&=&{\mathrm{exp}}(-x^2)(1+x^2+\frac{x^4}{2}).\end{aligned}
          :label: eqnljpmerealdistance

The above methodology works fine as long as the dispersion parameters
can be combined geometrically (:eq:`eqn. %s <eqncomb>`) in the same way as the
charges for electrostatics

.. math:: C^{ij}_{6,\mathrm{geom}} = \left(C^{ii}_6 \, C^{jj}_6\right)^{1/2}
          :label: eqnljpmegeom

For Lorentz-Berthelot combination rules (:eq:`eqn. %s <eqnlorentzberthelot>`),
the reciprocal part of this sum has to be calculated seven times due to
the splitting of the dispersion parameter according to

.. math:: C^{ij}_{6,\mathrm{L-B}} = (\sigma_i+\sigma_j)^6=\sum_{n=0}^{6} P_{n}\sigma_{i}^{n}\sigma_{j}^{(6-n)},
          :label: eqnljpmelorenztberthelot

for :math:`P_{n}` the Pascal triangle coefficients. This introduces a
non-negligible cost to the reciprocal part, requiring seven separate
FFTs, and therefore this has been the limiting factor in previous
attempts to implement LJ-PME. A solution to this problem is to use
geometrical combination rules in order to calculate an approximate
interaction parameter for the reciprocal part of the potential, yielding
a total interaction of

.. math:: \begin{aligned}
          V(r<r_c) & = & \underbrace{C^{\mathrm{dir}}_6 g(\beta r) r^{-6}}_{\mathrm{Direct \  space}} + \underbrace{C^\mathrm{recip}_{6,\mathrm{geom}} [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}} \nonumber \\
          &=& C^\mathrm{recip}_{6,\mathrm{geom}}r^{-6} + \left(C^{\mathrm{dir}}_6-C^\mathrm{recip}_{6,\mathrm{geom}}\right)g(\beta r)r^{-6} \\
          V(r>r_c) & = & \underbrace{C^\mathrm{recip}_{6,\mathrm{geom}} [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}}.\end{aligned}
          :label: eqnpmearith

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

.. math:: \begin{aligned} 
          V(r>r_c) &=& C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}.\end{aligned}
          :label: eqnljpmecorr3

For the case when :math:`C^{\mathrm{dir}}_6 \neq C^\mathrm{recip}_6`
this will retain an unmodified LJ force up to the cut-off, and the error
is an order of magnitude smaller than in simulations where the
direct-space interactions do not account for the approximation used in
reciprocal space. When using a VdW interaction modifier of
potential-shift, the constant

.. math:: \left(-C^{\mathrm{dir}}_6 + C^\mathrm{recip}_6 [1 - g(\beta r_c)]\right) r_c^{-6}
          :label: eqnljpmeconstant

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
