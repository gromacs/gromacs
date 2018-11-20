Enforced Rotation
-----------------

This module can be used to enforce the rotation of a group of atoms, as
*e.g.* a protein subunit. There are a variety of rotation potentials,
among them complex ones that allow flexible adaptations of both the
rotated subunit as well as the local rotation axis during the
simulation. An example application can be found in ref.
:ref:`145 <refKutzner2011>`.

.. _fig-rotation:

.. figure:: plots/rotation.*
   :width: 13.00000cm

   Comparison of fixed and flexible axis rotation. A:
   Rotating the sketched shape inside the white tubular cavity can
   create artifacts when a fixed rotation axis (dashed) is used. More
   realistically, the shape would revolve like a flexible pipe-cleaner
   (dotted) inside the bearing (gray). B: Fixed rotation
   around an axis :math:`\mathbf{v}` with a pivot point
   specified by the vector :math:`\mathbf{u}`.
   C: Subdividing the rotating fragment into slabs with
   separate rotation axes (:math:`\uparrow`) and pivot points
   (:math:`\bullet`) for each slab allows for flexibility. The distance
   between two slabs with indices :math:`n` and :math:`n+1` is
   :math:`\Delta x`.

.. _fig-equipotential:

.. figure:: plots/equipotential.*
   :width: 13.00000cm

   Selection of different rotation potentials and definition of
   notation. All four potentials :math:`V` (color coded) are shown for a
   single atom at position :math:`\mathbf{x}_j(t)`.
   A: Isotropic potential :math:`V^\mathrm{iso}`,
   B: radial motion potential :math:`V^\mathrm{rm}` and
   flexible potential :math:`V^\mathrm{flex}`, C–D: radial
   motion2 potential :math:`V^\mathrm{rm2}` and flexible2 potential
   :math:`V^\mathrm{flex2}` for :math:`\epsilon'\mathrm{ = }0\mathrm{ nm}^2`
   (C) and :math:`\epsilon'\mathrm{ = }0.01\mathrm{nm}^2`
   (D). The rotation axis is perpendicular to the plane
   and marked by :math:`\otimes`. The light gray contours indicate
   Boltzmann factors :math:`e^{-V/(k_B T)}` in the
   :math:`\mathbf{x}_j`-plane for :math:`T=300`\ K and
   :math:`k\mathrm{ = }200\mathrm{kJ}/(\mathrm{mol }\cdot\mathrm{nm}^2)`. The green
   arrow shows the direction of the force
   :math:`\mathbf{F}_{\!j}` acting on atom :math:`j`; the
   blue dashed line indicates the motion of the reference position.

Fixed Axis Rotation
^^^^^^^^^^^^^^^^^^^

Stationary Axis with an Isotropic Potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the fixed axis approach (see :numref:`Fig. %s B <fig-rotation>`),
torque on a group of :math:`N` atoms with positions
:math:`\mathbf{x}_i` (denoted “rotation group”) is applied
by rotating a reference set of atomic positions – usually their initial
positions :math:`\mathbf{y}_i^0` – at a constant angular
velocity :math:`\omega` around an axis defined by a direction vector
:math:`\hat{\mathbf{v}}` and a pivot point
:math:`\mathbf{u}`. To that aim, each atom with
position :math:`\mathbf{x}_i` is attracted by a “virtual
spring” potential to its moving reference position
:math:`\mathbf{y}_i = \mathbf{\Omega}(t) (\mathbf{y}_i^0 - \mathbf{u})`,
where :math:`\mathbf{\Omega}(t)` is a matrix that describes the rotation
around the axis. In the simplest case, the “springs” are described by a
harmonic potential,

.. math:: V^\mathrm{iso} = \frac{k}{2} \sum_{i=1}^{N} w_i \left[ \mathbf{\Omega}(t)
          (\mathbf{y}_i^0 - \mathbf{u}) - (\mathbf{x}_i - \mathbf{u})  \right]^2
          :label: eqnpotiso

with optional mass-weighted prefactors :math:`w_i = N \, m_i/M` with
total mass :math:`M = \sum_{i=1}^N m_i`. The rotation matrix
:math:`\mathbf{\Omega}(t)` is

.. math:: \mathbf{\Omega}(t) =  
          \left(   
          \begin{array}{ccc}
          \cos\omega t + v_x^2{\,\xi\,}& v_x v_y{\,\xi\,}- v_z\sin\omega t  & v_x v_z{\,\xi\,}+ v_y\sin\omega t\\
          v_x v_y{\,\xi\,}+ v_z\sin\omega t  & \cos\omega t + v_y^2{\,\xi\,}& v_y v_z{\,\xi\,}- v_x\sin\omega t\\
          v_x v_z{\,\xi\,}- v_y\sin\omega t  & v_y v_z{\,\xi\,}+ v_x\sin\omega t  & \cos\omega t + v_z^2{\,\xi\,}\\
          \end{array}
          \right)
          :label: eqnrotmat

where :math:`v_x`, :math:`v_y`, and :math:`v_z` are the components of
the normalized rotation vector :math:`\hat{\mathbf{v}}`,
and :math:`{\,\xi\,}:= 1-\cos(\omega t)`. As illustrated in
:numref:`Fig.  %s A <fig-equipotential>` for a single atom :math:`j`,
the rotation matrix :math:`\mathbf{\Omega}(t)` operates on the initial
reference positions
:math:`\mathbf{y}_j^0 = \mathbf{x}_j(t_0)`
of atom :math:`j` at :math:`t=t_0`. At a later time :math:`t`, the
reference position has rotated away from its initial place (along the
blue dashed line), resulting in the force

.. math:: \mathbf{F}_{\!j}^\mathrm{iso} 
          = -\nabla_{\!j} \, V^\mathrm{iso} 
          = k \, w_j \left[
          \mathbf{\Omega}(t) (\mathbf{y}_j^0 - \mathbf{u}) - (\mathbf{x}_j - \mathbf{u} ) \right]
          :label: eqnforcefixed

which is directed towards the reference position.

Pivot-Free Isotropic Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of a fixed pivot vector :math:`\mathbf{u}` this
potential uses the center of mass :math:`\mathbf{x}_c` of
the rotation group as pivot for the rotation axis,

.. math:: \mathbf{x}_c   = \frac{1}{M} \sum_{i=1}^N m_i \mathbf{x}_i 
          \mbox{and}
          \mathbf{y}_c^0 = \frac{1}{M} \sum_{i=1}^N m_i \mathbf{y}_i^0 \ ,
          :label: eqncom

which yields the “pivot-free” isotropic potential

.. math:: V^\mathrm{iso-pf} = \frac{k}{2} \sum_{i=1}^{N} w_i \left[ \mathbf{\Omega}(t)
          (\mathbf{y}_i^0 - \mathbf{y}_c^0) - (\mathbf{x}_i - \mathbf{x}_c) \right]^2 ,
          :label: eqnpotisopf

with forces

.. math:: \mathbf{F}_{\!j}^\mathrm{iso-pf} = k \, w_j 
          \left[ 
          \mathbf{\Omega}(t) ( \mathbf{y}_j^0 - \mathbf{y}_c^0) 
                           - ( \mathbf{x}_j   - \mathbf{x}_c )
          \right] .
          :label: eqnforceisopf

Without mass-weighting, the pivot :math:`\mathbf{x}_c` is
the geometrical center of the group.

Parallel Motion Potential Variant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The forces generated by the isotropic potentials
(eqns. :eq:`%s <eqnpotiso>` and :eq:`%s <eqnpotisopf>`) also contain components parallel to the
rotation axis and thereby restrain motions along the axis of either the
whole rotation group (in case of :math:`V^\mathrm{iso}`) or within the
rotation group, in case of :math:`V^\mathrm{iso-pf}`.
        
For cases where unrestrained motion along the axis is preferred, we have implemented a
“parallel motion” variant by eliminating all components parallel to the
rotation axis for the potential. This is achieved by projecting the
distance vectors between reference and actual positions

.. math:: \mathbf{r}_i = \mathbf{\Omega}(t) (\mathbf{y}_i^0 - \mathbf{u}) - (\mathbf{x}_i - \mathbf{u})
          :label: eqnrotdistvectors

onto the plane perpendicular to the rotation vector,

.. math:: \mathbf{r}_i^\perp :=  \mathbf{r}_i - (\mathbf{r}_i \cdot \hat{\mathbf{v}})\hat{\mathbf{v}}
          :label: eqnproject

yielding

.. math:: \begin{aligned}
          \nonumber
          V^\mathrm{pm} &=& \frac{k}{2} \sum_{i=1}^{N} w_i ( \mathbf{r}_i^\perp )^2 \\
                  &=& \frac{k}{2} \sum_{i=1}^{N} w_i
           \left\lbrace
           \mathbf{\Omega}(t)
             (\mathbf{y}_i^0 - \mathbf{u}) - (\mathbf{x}_i - \mathbf{u})  \right. \nonumber \\
          && \left. - \left\lbrace
          \left[ \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{u}) - (\mathbf{x}_i - \mathbf{u}) \right] \cdot\hat{\mathbf{v}}
            \right\rbrace\hat{\mathbf{v}} \right\rbrace^2
          \end{aligned}
          :label: eqnpotpm

and similarly

.. math:: \mathbf{F}_{\!j}^\mathrm{pm} = k \, w_j \, \mathbf{r}_j^\perp
          :label: eqnforcepm

Pivot-Free Parallel Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Replacing in eqn. :eq:`%s <eqnpotpm>` the fixed pivot
:math:`\mathbf{u}` by the center of mass
:math:`\mathbf{x_c}` yields the pivot-free variant of the
parallel motion potential. With

.. math:: \mathbf{s}_i = \mathbf{\Omega}(t) (\mathbf{y}_i^0 - \mathbf{y}_c^0) - (\mathbf{x}_i - \mathbf{x}_c)
          :label: eqnparrallelpotential

the respective potential and forces are

.. math:: \begin{aligned}
          V^\mathrm{pm-pf} &=& \frac{k}{2} \sum_{i=1}^{N} w_i ( \mathbf{s}_i^\perp )^2 \end{aligned}
          :label: eqnpotpmpf

.. math:: \begin{aligned}
          \mathbf{F}_{\!j}^\mathrm{pm-pf} &=& k \, w_j \, \mathbf{s}_j^\perp
          \end{aligned}
          :label: eqnforcepmpf

Radial Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^

In the above variants, the minimum of the rotation potential is either a
single point at the reference position
:math:`\mathbf{y}_i` (for the isotropic potentials) or a
single line through :math:`\mathbf{y}_i` parallel to the
rotation axis (for the parallel motion potentials). As a result, radial
forces restrict radial motions of the atoms. The two subsequent types of
rotation potentials, :math:`V^\mathrm{rm}` and :math:`V^\mathrm{rm2}`, drastically
reduce or even eliminate this effect. The first variant, :math:`V^\mathrm{rm}`
(:numref:`Fig. %s B <fig-equipotential>`), eliminates all force
components parallel to the vector connecting the reference atom and the
rotation axis,

.. math:: V^\mathrm{rm} = \frac{k}{2} \sum_{i=1}^{N} w_i \left[
          \mathbf{p}_i
          \cdot(\mathbf{x}_i - \mathbf{u}) \right]^2 ,
          :label: eqnpotrm

with

.. math::   \mathbf{p}_i := 
            \frac{\hat{\mathbf{v}}\times \mathbf{\Omega}(t) (\mathbf{y}_i^0 - \mathbf{u})} {\| \hat{\mathbf{v}}\times \mathbf{\Omega}(t) (\mathbf{y}_i^0 - \mathbf{u})\|} \ .
            :label: eqnpotrmpart2

This variant depends only on the distance
:math:`\mathbf{p}_i \cdot (\mathbf{x}_i -
\mathbf{u})` of atom :math:`i` from the plane spanned by
:math:`\hat{\mathbf{v}}` and
:math:`\mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{u})`.
The resulting force is

.. math:: \mathbf{F}_{\!j}^\mathrm{rm} =
           -k \, w_j \left[ \mathbf{p}_j\cdot(\mathbf{x}_j - \mathbf{u}) \right] \,\mathbf{p}_j \,  .
          :label: eqnpotrmforce

Pivot-Free Radial Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Proceeding similar to the pivot-free isotropic potential yields a
pivot-free version of the above potential. With

.. math:: \mathbf{q}_i := 
          \frac{\hat{\mathbf{v}}\times \mathbf{\Omega}(t) (\mathbf{y}_i^0 - \mathbf{y}_c^0)} {\| \hat{\mathbf{v}}\times \mathbf{\Omega}(t) (\mathbf{y}_i^0 - \mathbf{y}_c^0)\|} \, ,
          :label: eqnpotrmpfpart1

the potential and force for the pivot-free variant of the radial motion
potential read

.. math:: \begin{aligned}
          V^\mathrm{rm-pf} & = & \frac{k}{2} \sum_{i=1}^{N} w_i \left[
          \mathbf{q}_i
          \cdot(\mathbf{x}_i - \mathbf{x}_c)
          \right]^2 \, , \end{aligned}
          :label: eqnpotrmpf

.. math:: \begin{aligned}       
          \mathbf{F}_{\!j}^\mathrm{rm-pf} & = &
           -k \, w_j \left[ \mathbf{q}_j\cdot(\mathbf{x}_j - \mathbf{x}_c) \right] \,\mathbf{q}_j 
           + k   \frac{m_j}{M} \sum_{i=1}^{N} w_i \left[
           \mathbf{q}_i\cdot(\mathbf{x}_i - \mathbf{x}_c) \right]\,\mathbf{q}_i \, .
          \end{aligned}
          :label: eqnpotrmpfforce

Radial Motion 2 Alternative Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As seen in :numref:`Fig. %s B <fig-equipotential>`, the force
resulting from :math:`V^\mathrm{rm}` still contains a small, second-order
radial component. In most cases, this perturbation is tolerable; if not,
the following alternative, :math:`V^\mathrm{rm2}`, fully eliminates the
radial contribution to the force, as depicted in
:numref:`Fig. %s C <fig-equipotential>`,

.. math:: V^\mathrm{rm2} = 
          \frac{k}{2} \sum_{i=1}^{N} w_i\, 
          \frac{\left[ (\hat{\mathbf{v}} \times ( \mathbf{x}_i - \mathbf{u} ))
          \cdot \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{u}) \right]^2}
          {\| \hat{\mathbf{v}} \times (\mathbf{x}_i - \mathbf{u}) \|^2 +
          \epsilon'} \, ,
          :label: eqnpotrm2

where a small parameter :math:`\epsilon'` has been introduced to avoid
singularities. For :math:`\epsilon'\mathrm{ = }0\mathrm{nm}^2`, the
equipotential planes are spanned by :math:`\mathbf{x}_i -
\mathbf{u}` and :math:`\hat{\mathbf{v}}`,
yielding a force perpendicular to
:math:`\mathbf{x}_i - \mathbf{u}`, thus not
contracting or expanding structural parts that moved away from or toward
the rotation axis.

Choosing a small positive :math:`\epsilon'` (*e.g.*,
:math:`\epsilon'\mathrm{ = }0.01\mathrm{nm}^2`,
:numref:`Fig. %s D <fig-equipotential>`) in the denominator of
eqn. :eq:`%s <eqnpotrm2>` yields a well-defined potential and
continuous forces also close to the rotation axis, which is not the case
for :math:`\epsilon'\mathrm{ = }0\mathrm{nm}^2`
(:numref:`Fig. %s C <fig-equipotential>`). With

.. math:: \begin{aligned}
          \mathbf{r}_i & := & \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{u})\\
          \mathbf{s}_i & := & \frac{\hat{\mathbf{v}} \times (\mathbf{x}_i -
          \mathbf{u} ) }{ \| \hat{\mathbf{v}} \times (\mathbf{x}_i - \mathbf{u})
          \| } \equiv \; \Psi_{i} \;\; {\hat{\mathbf{v}} \times
          (\mathbf{x}_i-\mathbf{u} ) }\\
          \Psi_i^{*}   & := & \frac{1}{ \| \hat{\mathbf{v}} \times
          (\mathbf{x}_i-\mathbf{u}) \|^2 + \epsilon'}\end{aligned}
          :label: eqnpotrm2forcepart1

the force on atom :math:`j` reads

.. math:: \mathbf{F}_{\!j}^\mathrm{rm2}  = 
          - k\; 
          \left\lbrace w_j\;
          (\mathbf{s}_j\cdot\mathbf{r}_{\!j})\;
          \left[ \frac{\Psi_{\!j}^*   }{\Psi_{\!j}  }  \mathbf{r}_{\!j} 
               - \frac{\Psi_{\!j}^{ * 2}}{\Psi_{\!j}^3}
               (\mathbf{s}_j\cdot\mathbf{r}_{\!j})\mathbf{s}_j \right]
          \right\rbrace \times \hat{\mathbf{v}} .
          :label: eqnpotrm2force

Pivot-Free Radial Motion 2 Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pivot-free variant of the above potential is

.. math:: V{^\mathrm{rm2-pf}}= 
          \frac{k}{2} \sum_{i=1}^{N} w_i\, 
          \frac{\left[ (\hat{\mathbf{v}} \times ( \mathbf{x}_i - \mathbf{x}_c ))
          \cdot \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{y}_c) \right]^2}
          {\| \hat{\mathbf{v}} \times (\mathbf{x}_i - \mathbf{x}_c) \|^2 +
          \epsilon'} \, .
          :label: eqnpotrm2pf

With

.. math:: \begin{aligned}
          \mathbf{r}_i & := & \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{y}_c)\\
          \mathbf{s}_i & := & \frac{\hat{\mathbf{v}} \times (\mathbf{x}_i -
          \mathbf{x}_c ) }{ \| \hat{\mathbf{v}} \times (\mathbf{x}_i - \mathbf{x}_c)
          \| } \equiv \; \Psi_{i} \;\; {\hat{\mathbf{v}} \times
          (\mathbf{x}_i-\mathbf{x}_c ) }\\ \Psi_i^{*}   & := & \frac{1}{ \| \hat{\mathbf{v}} \times
          (\mathbf{x}_i-\mathbf{x}_c) \|^2 + \epsilon'}\end{aligned}
          :label: eqnpotrm2pfpart2

the force on atom :math:`j` reads

.. math:: \begin{aligned}
          \nonumber
          \mathbf{F}_{\!j}{^\mathrm{rm2-pf}}& = &
          - k\; 
          \left\lbrace w_j\;
          (\mathbf{s}_j\cdot\mathbf{r}_{\!j})\;
          \left[ \frac{\Psi_{\!j}^*   }{\Psi_{\!j}  } \mathbf{r}_{\!j} 
               - \frac{\Psi_{\!j}^{ * 2}}{\Psi_{\!j}^3}
               (\mathbf{s}_j\cdot\mathbf{r}_{\!j})\mathbf{s}_j \right]
          \right\rbrace \times \hat{\mathbf{v}}\\
               & &
          + k\;\frac{m_j}{M} \left\lbrace \sum_{i=1}^{N}
          w_i\;(\mathbf{s}_i\cdot\mathbf{r}_i) \; 
          \left[ \frac{\Psi_i^*   }{\Psi_i  }  \mathbf{r}_i
               - \frac{\Psi_i^{ * 2}}{\Psi_i^3} (\mathbf{s}_i\cdot\mathbf{r}_i )\;
               \mathbf{s}_i \right] \right\rbrace \times \hat{\mathbf{v}} \, .
          \end{aligned}
          :label: eqnpotrm2pfforce

Flexible Axis Rotation
~~~~~~~~~~~~~~~~~~~~~~

As sketched in :numref:`Fig. %s <fig-rotation>` A–B, the rigid body
behavior of the fixed axis rotation scheme is a drawback for many
applications. In particular, deformations of the rotation group are
suppressed when the equilibrium atom positions directly depend on the
reference positions. To avoid this limitation,
eqns. :eq:`%s <eqnpotrmpf>` and :eq:`%s <eqnpotrm2pf>`
will now be generalized towards a “flexible axis” as sketched in
:numref:`Fig. %s C <fig-rotation>`. This will be achieved by
subdividing the rotation group into a set of equidistant slabs
perpendicular to the rotation vector, and by applying a separate
rotation potential to each of these slabs.
:numref:`Fig. %s C <fig-rotation>` shows the midplanes of the slabs
as dotted straight lines and the centers as thick black dots.

To avoid discontinuities in the potential and in the forces, we define
“soft slabs” by weighing the contributions of each slab :math:`n` to the
total potential function :math:`V^\mathrm{flex}` by a Gaussian function

.. math:: g_n(\mathbf{x}_i) = \Gamma \ \mbox{exp} \left(
          -\frac{\beta_n^2(\mathbf{x}_i)}{2\sigma^2}  \right) ,
          :label: eqngaussian

centered at the midplane of the :math:`n`\ th slab. Here :math:`\sigma`
is the width of the Gaussian function, :math:`\Delta x` the distance
between adjacent slabs, and

.. math:: \beta_n(\mathbf{x}_i) := \mathbf{x}_i \cdot \hat{\mathbf{v}} - n \, \Delta x \, .
          :label: eqngaussianpart2

.. _fig-gaussian:

.. figure:: plots/gaussians.*
   :width: 6.50000cm

   Gaussian functions :math:`g_n` centered at :math:`n \, \Delta x` for
   a slab distance :math:`\Delta x = 1.5` nm and :math:`n \geq -2`.
   Gaussian function :math:`g_0` is highlighted in bold; the dashed line
   depicts the sum of the shown Gaussian functions.

A most convenient choice is :math:`\sigma = 0.7 \Delta x` and

.. math:: 1/\Gamma = \sum_{n \in Z}
          \mbox{exp}
          \left(-\frac{(n - \frac{1}{4})^2}{2\cdot 0.7^2}\right)
          \approx 1.75464 \, ,
          :label: eqngaussianpart3

which yields a nearly constant sum, essentially independent of
:math:`\mathbf{x}_i` (dashed line in
:numref:`Fig. %s <fig-gaussian>`), *i.e.*,

.. math:: \sum_{n \in Z} g_n(\mathbf{x}_i) =  1 + \epsilon(\mathbf{x}_i) \, ,
          :label: eqnnormal

with
:math:`| \epsilon(\mathbf{x}_i) | < 1.3\cdot 10^{-4}`.
This choice also implies that the individual contributions to the force
from the slabs add up to unity such that no further normalization is
required.

To each slab center :math:`\mathbf{x}_c^n`, all atoms
contribute by their Gaussian-weighted (optionally also mass-weighted)
position vectors
:math:`g_n(\mathbf{x}_i) \, \mathbf{x}_i`.
The instantaneous slab centers :math:`\mathbf{x}_c^n` are
calculated from the current positions
:math:`\mathbf{x}_i`,

.. math::  \mathbf{x}_c^n =
           \frac{\sum_{i=1}^N g_n(\mathbf{x}_i) \, m_i \, \mathbf{x}_i}
                {\sum_{i=1}^N g_n(\mathbf{x}_i) \, m_i} \, ,\\
           :label: eqndefx0 

while the reference centers :math:`\mathbf{y}_c^n` are
calculated from the reference positions
:math:`\mathbf{y}_i^0`,

.. math:: \mathbf{y}_c^n =
          \frac{\sum_{i=1}^N g_n(\mathbf{y}_i^0) \, m_i \, \mathbf{y}_i^0}
               {\sum_{i=1}^N g_n(\mathbf{y}_i^0) \, m_i} \, .
          :label: eqndefy0

Due to the rapid decay of :math:`g_n`, each slab will essentially
involve contributions from atoms located within :math:`\approx
3\Delta x` from the slab center only.

Flexible Axis Potential
^^^^^^^^^^^^^^^^^^^^^^^

We consider two flexible axis variants. For the first variant, the slab
segmentation procedure with Gaussian weighting is applied to the radial
motion potential
(eqn. :eq:`%s <eqnpotrmpf>` / :numref:`Fig. %s B <fig-equipotential>`),
yielding as the contribution of slab :math:`n`

.. math::  V^n = 
           \frac{k}{2} \sum_{i=1}^{N} w_i \, g_n(\mathbf{x}_i) 
           \left[
           \mathbf{q}_i^n
           \cdot
            (\mathbf{x}_i - \mathbf{x}_c^n) 
           \right]^2  ,
           :label: eqnflexpot

and a total potential function

.. math:: V^\mathrm{flex} = \sum_n V^n \, .
          :label: eqnpotflex

Note that the global center of mass :math:`\mathbf{x}_c`
used in eqn. :eq:`%s <eqnpotrmpf>` is now replaced by
:math:`\mathbf{x}_c^n`, the center of mass of the slab.
With

.. math:: \begin{aligned}
          \mathbf{q}_i^n & := & \frac{\hat{\mathbf{v}} \times
          \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{y}_c^n) }{ \| \hat{\mathbf{v}}
          \times \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{y}_c^n) \| } \\
          b_i^n         & := & \mathbf{q}_i^n \cdot (\mathbf{x}_i - \mathbf{x}_c^n) \, ,\end{aligned}
          :label: eqnflexpotpart2

the resulting force on atom :math:`j` reads

.. math:: \begin{aligned}
          \nonumber\hspace{-15mm}
          \mathbf{F}_{\!j}^\mathrm{flex} &=&
          - \, k \, w_j \sum_n g_n(\mathbf{x}_j) \, b_j^n \left\lbrace  \mathbf{q}_j^n -
          b_j^n \frac{\beta_n(\mathbf{x}_j)}{2\sigma^2} \hat{\mathbf{v}} \right\rbrace \\ & &
          + \, k \, m_j \sum_n \frac{g_n(\mathbf{x}_j)}{\sum_h g_n(\mathbf{x}_h)}
          \sum_{i=1}^{N} w_i \, g_n(\mathbf{x}_i) \, b_i^n \left\lbrace 
          \mathbf{q}_i^n -\frac{\beta_n(\mathbf{x}_j)}{\sigma^2}
          \left[ \mathbf{q}_i^n \cdot (\mathbf{x}_j - \mathbf{x}_c^n )\right]
          \hat{\mathbf{v}} \right\rbrace .
          \end{aligned}
          :label: eqnpotflexforce

Note that for :math:`V^\mathrm{flex}`, as defined, the slabs are fixed in
space and so are the reference centers
:math:`\mathbf{y}_c^n`. If during the simulation the
rotation group moves too far in :math:`\mathbf{v}`
direction, it may enter a region where – due to the lack of nearby
reference positions – no reference slab centers are defined, rendering
the potential evaluation impossible. We therefore have included a
slightly modified version of this potential that avoids this problem by
attaching the midplane of slab :math:`n=0` to the center of mass of the
rotation group, yielding slabs that move with the rotation group. This
is achieved by subtracting the center of mass
:math:`\mathbf{x}_c` of the group from the positions,

.. math:: \tilde{\mathbf{x}}_i = \mathbf{x}_i - \mathbf{x}_c \, , \mbox{\ \ \ and \ \ } 
          \tilde{\mathbf{y}}_i^0 = \mathbf{y}_i^0 - \mathbf{y}_c^0 \, ,
          :label: eqntrafo

such that

.. math:: \begin{aligned}
          V^\mathrm{flex-t} 
            & = & \frac{k}{2} \sum_n \sum_{i=1}^{N} w_i \, g_n(\tilde{\mathbf{x}}_i)
            \left[ \frac{\hat{\mathbf{v}} \times \mathbf{\Omega}(t)(\tilde{\mathbf{y}}_i^0
            - \tilde{\mathbf{y}}_c^n) }{ \| \hat{\mathbf{v}} \times
          \mathbf{\Omega}(t)(\tilde{\mathbf{y}}_i^0 -
          \tilde{\mathbf{y}}_c^n) \| }
          \cdot
           (\tilde{\mathbf{x}}_i - \tilde{\mathbf{x}}_c^n) 
          \right]^2 .
          \end{aligned}
          :label: eqnpotflext

To simplify the force derivation, and for efficiency reasons, we here
assume :math:`\mathbf{x}_c` to be constant, and thus
:math:`\partial \mathbf{x}_c / \partial x =
\partial \mathbf{x}_c / \partial y = \partial \mathbf{x}_c / \partial z = 0`.
The resulting force error is small (of order :math:`O(1/N)` or
:math:`O(m_j/M)` if mass-weighting is applied) and can therefore be
tolerated. With this assumption, the forces :math:`\mathbf{F}^\mathrm{flex-t}`
have the same form as eqn. :eq:`%s <eqnpotflexforce>`.

Flexible Axis 2 Alternative Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this second variant, slab segmentation is applied to
:math:`V^\mathrm{rm2}` (eqn. :eq:`%s <eqnpotrm2pf>`), resulting in
a flexible axis potential without radial force contributions
(:numref:`Fig. %s C <fig-equipotential>`),

.. math::   V{^\mathrm{flex2}}= 
            \frac{k}{2} \sum_{i=1}^{N} \sum_n w_i\,g_n(\mathbf{x}_i) 
            \frac{\left[ (\hat{\mathbf{v}} \times ( \mathbf{x}_i - \mathbf{x}_c^n ))
            \cdot \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{y}_c^n) \right]^2}
            {\| \hat{\mathbf{v}} \times (\mathbf{x}_i - \mathbf{x}_c^n) \|^2 +
            \epsilon'}
            :label: eqnpotflex2

With

.. math:: \begin{aligned}
          \mathbf{r}_i^n & := & \mathbf{\Omega}(t)(\mathbf{y}_i^0 - \mathbf{y}_c^n)\\
          \mathbf{s}_i^n & := & \frac{\hat{\mathbf{v}} \times (\mathbf{x}_i -
          \mathbf{x}_c^n ) }{ \| \hat{\mathbf{v}} \times (\mathbf{x}_i - \mathbf{x}_c^n)
          \| } \equiv \; \psi_{i} \;\; {\hat{\mathbf{v}} \times (\mathbf{x}_i-\mathbf{x}_c^n ) }\\
          \psi_i^{*}     & := & \frac{1}{ \| \hat{\mathbf{v}} \times (\mathbf{x}_i-\mathbf{x}_c^n) \|^2 + \epsilon'}\\
          W_j^n          & := & \frac{g_n(\mathbf{x}_j)\,m_j}{\sum_h g_n(\mathbf{x}_h)\,m_h}\\
          \mathbf{S}^n   & := & 
          \sum_{i=1}^{N} w_i\;g_n(\mathbf{x}_i)
          \; (\mathbf{s}_i^n\cdot\mathbf{r}_i^n)
          \left[ \frac{\psi_i^*   }{\psi_i  }  \mathbf{r}_i^n
               - \frac{\psi_i^{ * 2}}{\psi_i^3} (\mathbf{s}_i^n\cdot\mathbf{r}_i^n )\;
               \mathbf{s}_i^n \right] 
          \end{aligned}
          :label: eqnSn

the force on atom :math:`j` reads

.. math:: \begin{aligned}
          \nonumber
          \mathbf{F}_{\!j}{^\mathrm{flex2}}& = &
          - k\; 
          \left\lbrace \sum_n w_j\;g_n(\mathbf{x}_j)\;
          (\mathbf{s}_j^n\cdot\mathbf{r}_{\!j}^n)\;
          \left[ \frac{\psi_j^*   }{\psi_j  }  \mathbf{r}_{\!j}^n 
               - \frac{\psi_j^{ * 2}}{\psi_j^3} (\mathbf{s}_j^n\cdot\mathbf{r}_{\!j}^n)\;
               \mathbf{s}_{\!j}^n \right] \right\rbrace \times \hat{\mathbf{v}} \\
          \nonumber
          & &
          + k \left\lbrace \sum_n W_{\!j}^n \, \mathbf{S}^n \right\rbrace \times
          \hat{\mathbf{v}}
          - k \left\lbrace \sum_n W_{\!j}^n \; \frac{\beta_n(\mathbf{x}_j)}{\sigma^2} \frac{1}{\psi_j}\;\; 
          \mathbf{s}_j^n \cdot 
          \mathbf{S}^n \right\rbrace \hat{\mathbf{v}}\\ 
          & & 
          + \frac{k}{2} \left\lbrace \sum_n w_j\;g_n(\mathbf{x}_j)
          \frac{\beta_n(\mathbf{x}_j)}{\sigma^2} 
          \frac{\psi_j^*}{\psi_j^2}( \mathbf{s}_j^n \cdot \mathbf{r}_{\!j}^n )^2 \right\rbrace
          \hat{\mathbf{v}} .
          \end{aligned}
          :label: eqnpotflex2force

Applying transformation :eq:`%s <eqntrafo>` yields a
“translation-tolerant” version of the flexible2 potential,
:math:`V{^\mathrm{flex2 - t}}`. Again, assuming that
:math:`\partial \mathbf{x}_c / \partial x`,
:math:`\partial \mathbf{x}_c /
\partial y`, :math:`\partial \mathbf{x}_c / \partial z`
are small, the resulting equations for :math:`V{^\mathrm{flex2 - t}}`
and :math:`\mathbf{F}{^\mathrm{flex2 - t}}` are similar
to those of :math:`V^\mathrm{flex2}` and
:math:`\mathbf{F}^\mathrm{flex2}`.

Usage
~~~~~

To apply enforced rotation, the particles :math:`i` that are to be
subjected to one of the rotation potentials are defined via index groups
``rot-group0``, ``rot-group1``, etc., in the
:ref:`mdp` input file. The reference positions
:math:`\mathbf{y}_i^0` are read from a special
:ref:`trr` file provided to :ref:`grompp <gmx grompp>`. If no such
file is found, :math:`\mathbf{x}_i(t=0)` are used as
reference positions and written to :ref:`trr` such that they
can be used for subsequent setups. All parameters of the potentials such
as :math:`k`, :math:`\epsilon'`, etc.
(:numref:`Table %s <tab-vars>`) are provided as :ref:`mdp`
parameters; ``rot-type`` selects the type of the potential.
The option ``rot-massw`` allows to choose whether or not to
use mass-weighted averaging. For the flexible potentials, a cutoff value
:math:`g_n^\mathrm{min}` (typically :math:`g_n^\mathrm{min}=0.001`)
makes sure that only significant contributions to :math:`V` and
:math:`\mathbf{F}` are evaluated, *i.e.* terms with
:math:`g_n(\mathbf{x}) < g_n^\mathrm{min}` are omitted.
:numref:`Table %s <tab-quantities>` summarizes observables that are
written to additional output files and which are described below.

.. |ROTISO| replace:: V\ :math:`^\mathrm{iso}`
.. |ROTISOPF| replace:: V\ :math:`^\mathrm{iso-pf}`
.. |ROTPM| replace:: V\ :math:`^\mathrm{pm}`
.. |ROTPMPF| replace:: V\ :math:`^\mathrm{pm-pf}`
.. |ROTRM| replace:: V\ :math:`^\mathrm{rm}`
.. |ROTRMPF| replace:: V\ :math:`^\mathrm{rm-pf}`
.. |ROTRM2| replace:: V\ :math:`^\mathrm{rm2}`
.. |ROTRM2PF| replace:: V\ :math:`^\mathrm{rm2-pf}`
.. |ROTFL| replace:: V\ :math:`^\mathrm{flex}`
.. |ROTFLT| replace:: V\ :math:`^\mathrm{flex-t}`
.. |ROTFL2| replace:: V\ :math:`^\mathrm{flex2}`
.. |ROTFLT2| replace:: V\ :math:`^\mathrm{flex2-t}`
.. |KUNIT| replace:: :math:`\frac{\mathrm{kJ}}{\mathrm{mol} \cdot \mathrm{nm}^2}`
.. |BFX| replace:: **x**
.. |KMA| replace:: :math:`k`
.. |VECV| replace:: :math:`\hat{\mathbf{v}}`
.. |VECU| replace:: :math:`\mathbf{u}`
.. |OMEG| replace:: :math:`\omega`
.. |EPSP| replace:: :math:`{\epsilon}'`
.. |DELX| replace:: :math:`{\Delta}x`
.. |GMIN| replace:: :math:`g_n^\mathrm{min}`
.. |CIPS| replace:: :math:`^\circ`\ /ps
.. |NM2| replace:: nm\ :math:`^2`
.. |REF1| replace:: \ :eq:`eqnpotiso`
.. |REF2| replace:: \ :eq:`eqnpotisopf`
.. |REF3| replace:: \ :eq:`eqnpotpm`
.. |REF4| replace:: \ :eq:`eqnpotpmpf`
.. |REF5| replace:: \ :eq:`eqnpotrm`
.. |REF6| replace:: \ :eq:`eqnpotrmpf`
.. |REF7| replace:: \ :eq:`eqnpotrm2`
.. |REF8| replace:: \ :eq:`eqnpotrm2pf`
.. |REF9| replace:: \ :eq:`eqnpotflex`
.. |REF10| replace:: \ :eq:`eqnpotflext`
.. |REF11| replace:: \ :eq:`eqnpotflex2`

.. _tab-vars:

.. table:: Parameters used by the various rotation potentials.
           |BFX| indicate which parameter is actually used for a given potential
           :widths: auto
           :align: center

           +------------------------------------------+---------+--------+--------+--------+--------+-----------+-----------+
           | parameter                                | |KMA|   | |VECV| | |VECU| | |OMEG| | |EPSP| | |DELX|    | |GMIN|    |
           +------------------------------------------+---------+--------+--------+--------+--------+-----------+-----------+
           | :ref:`mdp` input variable name           | k       | vec    | pivot  | rate   | eps    | slab-dist | min-gauss |
           +------------------------------------------+---------+--------+--------+--------+--------+-----------+-----------+
           | unit                                     | |KUNIT| | ``-``  | nm     | |CIPS| | |NM2|  | nm        | ``-``     |
           +================================+=========+=========+========+========+========+========+===========+===========+
           | fixed axis potentials:         | eqn.                                                                          |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | isotropic         | |ROTISO|   | |REF1|  | |BFX|   | |BFX|  | |BFX|  | |BFX|  | ``-``  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | --- pivot-free    | |ROTISOPF| | |REF2|  | |BFX|   | |BFX|  | ``-``  | |BFX|  | ``-``  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | parallel motion   | |ROTPM|    | |REF3|  | |BFX|   | |BFX|  | |BFX|  | |BFX|  | ``-``  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | --- pivot-free    | |ROTPMPF|  | |REF4|  | |BFX|   | |BFX|  | ``-``  | |BFX|  | ``-``  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | radial motion     | |ROTRM|    | |REF5|  | |BFX|   | |BFX|  | |BFX|  | |BFX|  | ``-``  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | --- pivot-free    | |ROTRMPF|  | |REF6|  | |BFX|   | |BFX|  | ``-``  | |BFX|  | ``-``  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | radial motion 2   | |ROTRM2|   | |REF7|  | |BFX|   | |BFX|  | |BFX|  | |BFX|  | |BFX|  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | --- pivot-free    | |ROTRM2PF| | |REF8|  | |BFX|   | |BFX|  | ``-``  | |BFX|  | |BFX|  | ``-``     | ``-``     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | flexible axis potentials:      | eqn.                                                                          | 
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | flexible          | |ROTFL|    | |REF9|  | |BFX|   | |BFX|  | ``-``  | |BFX|  | ``-``  | |BFX|     | |BFX|     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | --- transl. tol   | |ROTFLT|   | |REF10| | |BFX|   | |BFX|  | ``-``  | |BFX|  | ``-``  | |BFX|     | |BFX|     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | flexible 2        | |ROTFL2|   | |REF11| | |BFX|   | |BFX|  | ``-``  | |BFX|  | |BFX|  | |BFX|     | |BFX|     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+
           | --- transl. tol   | |ROTFLT2|  | ``-``   | |BFX|   | |BFX|  | ``-``  | |BFX|  | |BFX|  | |BFX|     | |BFX|     |
           +-------------------+------------+---------+---------+--------+--------+--------+--------+-----------+-----------+


| 

.. |VT|      replace:: :math:`V(t)`
.. |THET|    replace:: :math:`\theta_{	\mathrm{ref}}(t)`
.. |THETAV|  replace:: :math:`\theta_{\mathrm{av}}(t)`
.. |THETFIT| replace:: :math:`\theta_{\mathrm{fit}}(t)`, :math:`\theta_{\mathrm{fit}}(t,n)`
.. |YVEC|    replace:: :math:`\mathbf{y}_{0}(n)`, :math:`\mathbf{x}_{0}(t,n)`
.. |TAUT|    replace:: :math:`\tau(t)`
.. |TAUTN|   replace:: :math:`\tau(t,n)`
.. |REFT|  replace:: :numref:`see Table %s <tab-vars>`
.. |REFEQ| replace:: :math:`\theta_{\mathrm{ref}}(t)=\omega t`
.. |REF12| replace:: \ :eq:`eqnavangle`
.. |REF13| replace:: \ :eq:`eqnrmsdfit`
.. |REF14| replace:: \ :eq:`eqndefx0`\ ,\ :eq:`eqndefy0`
.. |REF15| replace:: \ :eq:`eqntorque` 

.. _tab-quantities:

.. table:: Quantities recorded in output files during enforced rotation simulations.
           All slab-wise data is written every ``nstsout`` steps, other rotation data every ``nstrout`` steps.
           :widths: auto
           :align: center

           +------------+---------+------------+--------------------+-------+----------+
           | quantity   | unit    | equation   | output file        | fixed | flexible |
           +============+=========+============+====================+=======+==========+
           | |VT|       | kJ/mol  | |REFT|     | ``rotation``       | |BFX| | |BFX|    |
           +------------+---------+------------+--------------------+-------+----------+
           | |THET|     | degrees | |REFEQ|    | ``rotation``       | |BFX| | |BFX|    |
           +------------+---------+------------+--------------------+-------+----------+
           | |THETAV|   | degrees | |REF12|    | ``rotation``       | |BFX| | ``-``    |
           +------------+---------+------------+--------------------+-------+----------+
           | |THETFIT|  | degrees | |REF13|    | ``rotangles``      | ``-`` | |BFX|    |
           +------------+---------+------------+--------------------+-------+----------+
           | |YVEC|     | nm      | |REF14|    | ``rotslabs``       | ``-`` | |BFX|    |
           +------------+---------+------------+--------------------+-------+----------+
           | |TAUT|     | kJ/mol  | |REF15|    | ``rotation``       | |BFX| | ``-``    |
           +------------+---------+------------+--------------------+-------+----------+
           | |TAUTN|    | kJ/mol  | |REF15|    | ``rottorque``      | ``-`` | |BFX|    |
           +------------+---------+------------+--------------------+-------+----------+




Angle of Rotation Groups: Fixed Axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For fixed axis rotation, the average angle :math:`\theta_\mathrm{av}(t)`
of the group relative to the reference group is determined via the
distance-weighted angular deviation of all rotation group atoms from
their reference positions,

.. math::  \theta_\mathrm{av} = \left. \sum_{i=1}^{N} r_i \ \theta_i \right/ \sum_{i=1}^N r_i \ .
           :label: eqnavangle

Here, :math:`r_i` is the distance of the reference position to the
rotation axis, and the difference angles :math:`\theta_i` are determined
from the atomic positions, projected onto a plane perpendicular to the
rotation axis through pivot point :math:`\mathbf{u}` (see
eqn. :eq:`%s <eqnproject>` for the definition of
:math:`\perp`),

.. math:: \cos \theta_i = 
          \frac{(\mathbf{y}_i-\mathbf{u})^\perp \cdot (\mathbf{x}_i-\mathbf{u})^\perp}
               { \| (\mathbf{y}_i-\mathbf{u})^\perp \cdot (\mathbf{x}_i-\mathbf{u})^\perp
               \| } \ .
          :label: eqnavanglepart2

The sign of :math:`\theta_\mathrm{av}` is chosen such that
:math:`\theta_\mathrm{av} > 0` if the actual structure rotates ahead of
the reference.

Angle of Rotation Groups: Flexible Axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For flexible axis rotation, two outputs are provided, the angle of the
entire rotation group, and separate angles for the segments in the
slabs. The angle of the entire rotation group is determined by an RMSD
fit of :math:`\mathbf{x}_i` to the reference positions
:math:`\mathbf{y}_i^0` at :math:`t=0`, yielding
:math:`\theta_\mathrm{fit}` as the angle by which the reference has to
be rotated around :math:`\hat{\mathbf{v}}` for the optimal
fit,

.. math::  \mathrm{RMSD} \big( \mathbf{x}_i,\ \mathbf{\Omega}(\theta_\mathrm{fit})
           \mathbf{y}_i^0 \big) \stackrel{!}{=} \mathrm{min} \, .
           :label: eqnrmsdfit

To determine the local angle for each slab :math:`n`, both reference
and actual positions are weighted with the Gaussian function of slab
:math:`n`, and :math:`\theta_\mathrm{fit}(t,n)` is calculated as in
eqn. :eq:`%s <eqnrmsdfit>` from the Gaussian-weighted
positions.

For all angles, the :ref:`mdp` input option
``rot-fit-method`` controls whether a normal RMSD fit is
performed or whether for the fit each position
:math:`\mathbf{x}_i` is put at the same distance to the
rotation axis as its reference counterpart
:math:`\mathbf{y}_i^0`. In the latter case, the RMSD
measures only angular differences, not radial ones.

Angle Determination by Searching the Energy Minimum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, for ``rot-fit-method = potential``, the angle
of the rotation group is determined as the angle for which the rotation
potential energy is minimal. Therefore, the used rotation potential is
additionally evaluated for a set of angles around the current reference
angle. In this case, the ``rotangles.log`` output file
contains the values of the rotation potential at the chosen set of
angles, while ``rotation.xvg`` lists the angle with minimal
potential energy.

Torque
^^^^^^

The torque :math:`\mathbf{\tau}(t)` exerted by the
rotation potential is calculated for fixed axis rotation via

.. math:: \mathbf{\tau}(t) = \sum_{i=1}^{N} \mathbf{r}_i(t) \times \mathbf{f}_{\!i}^\perp(t) ,
          :label: eqntorque

where :math:`\mathbf{r}_i(t)` is the distance vector from
the rotation axis to :math:`\mathbf{x}_i(t)` and
:math:`\mathbf{f}_{\!i}^\perp(t)` is the force component
perpendicular to :math:`\mathbf{r}_i(t)` and
:math:`\hat{\mathbf{v}}`. For flexible axis rotation,
torques :math:`\mathbf{\tau}_{\!n}` are calculated for
each slab using the local rotation axis of the slab and the
Gaussian-weighted positions.
