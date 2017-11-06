Special Topics
==============

Free energy implementation
--------------------------
For free energy calculations, there are two things that must be
specified; the end states, and the pathway connecting the end states.
The end states can be specified in two ways. The most straightforward is
through the specification of end states in the topology file. Most
potential forms support both an :math:`A` state and a :math:`B` state.
Whenever both states are specified, then the :math:`A` state corresponds
to the initial free energy state, and the :math:`B` state corresponds to
the final state.

In some cases, the end state can also be defined in some cases without
altering the topology, solely through the :ref:`mdp` file,
through the use of the
``couple-moltype``,
``couple-lambda0``,
``couple-lambda1``, and ``couple-intramol`` :ref:`mdp`
keywords. Any molecule type selected in ``couple-moltype``
will automatically have a :math:`B` state implicitly constructed (and
the :math:`A` state redefined) according to the
``couple-lambda`` keywords. ``couple-lambda0``
and ``couple-lambda1`` define the non-bonded parameters that
are present in the :math:`A` state (``couple-lambda0``) and
the :math:`B` state (``couple-lambda1``). The choices are
``q``,
``vdw``, and ``vdw-q``; these indicate the Coulombic, van der Waals, or
both parameters that are turned on in the respective state.

Once the end states are defined, then the path between the end states
has to be defined. This path is defined solely in the .mdp file.
Starting in 4.6, :math:`\lambda` is a vector of components, with
Coulombic, van der Waals, bonded, restraint, and mass components all
able to be adjusted independently. This makes it possible to turn off
the Coulombic term linearly, and then the van der Waals using soft core,
all in the same simulation. This is especially useful for replica
exchange or expanded ensemble simulations, where it is important to
sample all the way from interacting to non-interacting states in the
same simulation to improve sampling.

``fep-lambdas`` is the default array of :math:`\lambda`
values ranging from 0 to 1. All of the other lambda arrays use the
values in this array if they are not specified. The previous behavior,
where the pathway is controlled by a single :math:`\lambda` variable,
can be preserved by using only ``fep-lambdas`` to define the
pathway.

For example, if you wanted to first to change the Coulombic terms, then
the van der Waals terms, changing bonded at the same time rate as the
van der Waals, but changing the restraints throughout the first
two-thirds of the simulation, then you could use this :math:`\lambda`
vector:

::

    coul-lambdas           = 0.0 0.2 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    vdw-lambdas            = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
    bonded-lambdas         = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
    restraint-lambdas      = 0.0 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0

This is also equivalent to:

::

    fep-lambdas            = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
    coul-lambdas           = 0.0 0.2 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    restraint-lambdas      = 0.0 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0

The ``fep-lambda`` array, in this case, is being used as the
default to fill in the bonded and van der Waals :math:`\lambda` arrays.
Usually, it’s best to fill in all arrays explicitly, just to make sure
things are properly assigned.

If you want to turn on only restraints going from :math:`A` to
:math:`B`, then it would be:

::

    restraint-lambdas      = 0.0 0.1 0.2 0.4 0.6 1.0

and all of the other components of the :math:`\lambda` vector would be
left in the :math:`A` state.

To compute free energies with a vector :math:`\lambda` using
thermodynamic integration, then the TI equation becomes vector equation:

.. math:: \Delta F = \int \langle \nabla H \rangle \cdot d\vec{\lambda}

or for finite differences:

.. math:: \Delta F \approx \int \sum \langle \nabla H \rangle \cdot \Delta\lambda

The external `pymbar script <https://SimTK.org/home/pymbar>`__
can compute this integral automatically
from the |Gromacs| ``dhdl.xvg`` output.

Potential of mean force
-----------------------

A potential of mean force (PMF) is a potential that is obtained by
integrating the mean force from an ensemble of configurations. In
|Gromacs|, there are several different methods to calculate the mean
force. Each method has its limitations, which are listed below.

-  **pull code:** between the centers of mass of molecules or groups of
   molecules.

-  **free-energy code with harmonic bonds or constraints:** between
   single atoms.

-  **free-energy code with position restraints:** changing the
   conformation of a relatively immobile group of atoms.

-  **pull code in limited cases:** between groups of atoms that are part
   of a larger molecule for which the bonds are constrained with SHAKE
   or LINCS. If the pull group if relatively large, the pull code can be
   used.

The pull and free-energy code a described in more detail in the
following two sections.

Entropic effects
^^^^^^^^^^^^^^^^

When a distance between two atoms or the centers of mass of two groups
is constrained or restrained, there will be a purely entropic
contribution to the PMF due to the rotation of the two groups Neumann
(1980). For a system of two non-interacting masses the potential of mean
force is:

.. math:: V_{pmf}(r) = -(n_c - 1) k_B T \log(r)

where :math:`n_c` is the number of dimensions in which the constraint
works (i.e. :math:`n_c=3` for a normal constraint and :math:`n_c=1` when
only the :math:`z`-direction is constrained). Whether one needs to
correct for this contribution depends on what the PMF should represent.
When one wants to pull a substrate into a protein, this entropic term
indeed contributes to the work to get the substrate into the protein.
But when calculating a PMF between two solutes in a solvent, for the
purpose of simulating without solvent, the entropic contribution should
be removed. **Note** that this term can be significant; when at 300K the
distance is halved, the contribution is 3.5 kJ mol\ :math:`^{-1}`.

Non-equilibrium pulling
-----------------------

When the distance between two groups is changed continuously, work is
applied to the system, which means that the system is no longer in
equilibrium. Although in the limit of very slow pulling the system is
again in equilibrium, for many systems this limit is not reachable
within reasonable computational time. However, one can use the Jarzynski
relation Jarzynski () to obtain the equilibrium free-energy difference
:math:`\Delta G` between two distances from many non-equilibrium
simulations:

.. math:: \Delta G_{AB} = -k_BT \log \left\langle e^{-\beta W_{AB}} \right\rangle_A
          :label: eqJarz

where :math:`W_{AB}` is the work performed to force the system along
one path from state A to B, the angular bracket denotes averaging over a
canonical ensemble of the initial state A and :math:`\beta=1/k_B T`.

The pull code
-------------

The pull code applies forces or
constraints between the centers of mass of one or more pairs of groups
of atoms. Each pull reaction coordinate is called a “coordinate” and it
operates on usually two, but sometimes more, pull groups. A pull group
can be part of one or more pull coordinates. Furthermore, a coordinate
can also operate on a single group and an absolute reference position in
space. The distance between a pair of groups can be determined in 1, 2
or 3 dimensions, or can be along a user-defined vector. The reference
distance can be constant or can change linearly with time. Normally all
atoms are weighted by their mass, but an additional weighting factor can
also be used.

.. figure:: plots/pull.eps

   Schematic picture of pulling a lipid out of a lipid bilayer with
   umbrella pulling. :math:`V_{rup}` is the velocity at which the spring
   is retracted, :math:`Z_{link}` is the atom to which the spring is
   attached and :math:`Z_{spring}` is the location of the spring.

Several different pull types, i.e. ways to apply the pull force, are
supported, and in all cases the reference distance can be constant or
linearly changing with time.

#. **Umbrella pulling**
   A harmonic potential is applied between the centers of mass of two
   groups. Thus, the force is proportional to the displacement.

#. **Constraint pulling**
   The distance between the centers of mass of two groups is
   constrained. The constraint force can be written to a file. This
   method uses the SHAKE algorithm but only needs 1 iteration to be
   exact if only two groups are constrained.

#. **Constant force pulling** A constant force is applied between the
   centers of mass of two groups. Thus, the potential is linear. In this
   case there is no reference distance of pull rate.

#. **Flat bottom pulling** Like umbrella pulling, but the potential and
   force are zero for coordinate values below
   (``pull-coord?-type = flat-bottom``) or above
   (``pull-coord?-type = flat-bottom-high``) a reference
   value. This is useful for restraining e.g. the distance between two
   molecules to a certain region.

In addition, there are different types of reaction coordinates,
so-called pull geometries. These are set with the :ref:`mdp`
option ``pull-coord?-geometry``.

Definition of the center of mass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In |Gromacs|, there are three ways to define the center of mass of a
group. The standard way is a “plain” center of mass, possibly with
additional weighting factors. With periodic boundary conditions it is no
longer possible to uniquely define the center of mass of a group of
atoms. Therefore, a reference atom is used. For determining the center
of mass, for all other atoms in the group, the closest periodic image to
the reference atom is used. This uniquely defines the center of mass. By
default, the middle (determined by the order in the topology) atom is
used as a reference atom, but the user can also select any other atom if
it would be closer to center of the group.

For a layered system, for instance a lipid bilayer, it may be of
interest to calculate the PMF of a lipid as function of its distance
from the whole bilayer. The whole bilayer can be taken as reference
group in that case, but it might also be of interest to define the
reaction coordinate for the PMF more locally. The :ref:`mdp`
option ``pull-coord?-geometry = cylinder`` does not use all
the atoms of the reference group, but instead dynamically only those
within a cylinder with radius ``pull-cylinder-r`` around the
pull vector going through the pull group. This only works for distances
defined in one dimension, and the cylinder is oriented with its long
axis along this one dimension. To avoid jumps in the pull force,
contributions of atoms are weighted as a function of distance (in
addition to the mass weighting):

.. math::

   \begin{aligned}
   w(r < r_\mathrm{cyl}) & = &
   1-2 \left(\frac{r}{r_\mathrm{cyl}}\right)^2 + \left(\frac{r}{r_\mathrm{cyl}}\right)^4 \\
   w(r \geq r_\mathrm{cyl}) & = & 0\end{aligned}

Note that the radial dependence on the weight causes a radial force on
both cylinder group and the other pull group. This is an undesirable,
but unavoidable effect. To minimize this effect, the cylinder radius
should be chosen sufficiently large. The effective mass is 0.47 times
that of a cylinder with uniform weights and equal to the mass of uniform
cylinder of 0.79 times the radius.

.. figure:: plots/pullref.eps

   Comparison of a plain center of mass reference group versus a
   cylinder reference group applied to interface systems. C is the
   reference group. The circles represent the center of mass of two
   groups plus the reference group, :math:`d_c` is the reference
   distance.

For a group of molecules in a periodic system, a plain reference group
might not be well-defined. An example is a water slab that is connected
periodically in :math:`x` and :math:`y`, but has two liquid-vapor
interfaces along :math:`z`. In such a setup, water molecules can
evaporate from the liquid and they will move through the vapor, through
the periodic boundary, to the other interface. Such a system is
inherently periodic and there is no proper way of defining a “plain”
center of mass along :math:`z`. A proper solution is to using a cosine
shaped weighting profile for all atoms in the reference group. The
profile is a cosine with a single period in the unit cell. Its phase is
optimized to give the maximum sum of weights, including mass weighting.
This provides a unique and continuous reference position that is nearly
identical to the plain center of mass position in case all atoms are all
within a half of the unit-cell length. See ref O. Engin and Hess (2010)
for details.

When relative weights :math:`w_i` are used during the calculations,
either by supplying weights in the input or due to cylinder geometry or
due to cosine weighting, the weights need to be scaled to conserve
momentum:

.. math::

   w'_i = w_i
   \left. \sum_{j=1}^N w_j \, m_j \right/ \sum_{j=1}^N w_j^2 \, m_j

where :math:`m_j` is the mass of atom :math:`j` of the group. The mass
of the group, required for calculating the constraint force, is:

.. math:: M = \sum_{i=1}^N w'_i \, m_i

The definition of the weighted center of mass is:

.. math:: {\mbox{\boldmath ${r}$}}_{com} = \left. \sum_{i=1}^N w'_i \, m_i \, {\mbox{\boldmath ${r}$}}_i \right/ M

From the centers of mass the AFM, constraint, or umbrella force
:math:`{\mbox{\boldmath ${F}$}}_{\!com}` on each group can be
calculated. The force on the center of mass of a group is redistributed
to the atoms as follows:

.. math:: {\mbox{\boldmath ${F}$}}_{\!i} = \frac{w'_i \, m_i}{M} \, {\mbox{\boldmath ${F}$}}_{\!com}

Definition of the pull direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most common setup is to pull along the direction of the vector
containing the two pull groups, this is selected with
``pull-coord?-geometry = distance``. You might want to pull
along a certain vector instead, which is selected with
``pull-coord?-geometry = direction``. But this can cause
unwanted torque forces in the system, unless you pull against a
reference group with (nearly) fixed orientation, e.g. a membrane protein
embedded in a membrane along x/y while pulling along z. If your
reference group does not have a fixed orientation, you should probably
use ``pull-coord?-geometry = direction-relative``, see
:numref:`Fig. (%s) <fig-pulldirrel>`. Since the potential now depends
on the coordinates of two additional groups defining the orientation,
the torque forces will work on these two groups.

.. _fig-pulldirrel:

.. figure:: plots/pulldirrel.eps

   The pull setup for geometry ``direction-relative``. The
   “normal” pull groups are 1 and 2. Groups 3 and 4 define the pull
   direction and thus the direction of the normal pull forces (red).
   This leads to reaction forces (blue) on groups 3 and 4, which are
   perpendicular to the pull direction. Their magnitude is given by the
   “normal” pull force times the ratio of :math:`d_p` and the distance
   between groups 3 and 4.

Definition of the angle and dihedral pull geometries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Four pull groups are required for ``pull-coord?-geometry =
angle``. In the same way as for geometries with two groups, each
consecutive pair of groups :math:`i` and :math:`i+1` define a vector
connecting the COMs of groups :math:`i` and :math:`i+1`. The angle is
defined as the angle between the two resulting vectors. E.g., the
:ref:`mdp` option ``pull-coord?-groups = 1 2 2 4``
defines the angle between the vector from the COM of group 1 to the COM
of group 2 and the vector from the COM of group 2 to the COM of group 4.
The angle takes values in the closed interval [0, 180] deg. For
``pull-coord?-geometry = angle-axis`` the angle is defined
with respect to a reference axis given by
``pull-coord?-vec`` and only two groups need to be given.
The dihedral geometry requires six pull groups. These pair up in the
same way as described above and so define three vectors. The dihedral
angle is defined as the angle between the two planes spanned by the two
first and the two last vectors. Equivalently, the dihedral angle can be
seen as the angle between the first and the third vector when these
vectors are projected onto a plane normal to the second vector (the axis
vector). As an example, consider a dihedral angle involving four groups:
1, 5, 8 and 9. Here, the :ref:`mdp` option
``pull-coord?-groups = 8 1 1 5 5 9`` specifies the three
vectors that define the dihedral angle: the first vector is the COM
distance vector from group 8 to 1, the second vector is the COM distance
vector from group 1 to 5, and the third vector is the COM distance
vector from group 5 to 9. The dihedral angle takes values in the
interval (-180, 180] deg and has periodic boundaries.

Limitations
^^^^^^^^^^^

There is one theoretical limitation: strictly speaking, constraint
forces can only be calculated between groups that are not connected by
constraints to the rest of the system. If a group contains part of a
molecule of which the bond lengths are constrained, the pull constraint
and LINCS or SHAKE bond constraint algorithms should be iterated
simultaneously. This is not done in |Gromacs|. This means that for
simulations with ``constraints = all-bonds`` in the
:ref:`mdp` file pulling is, strictly speaking, limited to
whole molecules or groups of molecules. In some cases this limitation
can be avoided by using the free energy code, see
sec. 
. In practice, the errors caused by not
iterating the two constraint algorithms can be negligible when the pull
group consists of a large amount of atoms and/or the pull force is
small. In such cases, the constraint correction displacement of the pull
group is small compared to the bond lengths.

Enforced Rotation
-----------------

This module can be used to enforce the rotation of a group of atoms, as
*e.g.* a protein subunit. There are a variety of rotation potentials,
among them complex ones that allow flexible adaptations of both the
rotated subunit as well as the local rotation axis during the
simulation. An example application can be found in ref. Kutzner, Czub,
and Grubmüller (2011).

.. _fig-rotation:

.. figure:: plots/rotation.pdf

   Comparison of fixed and flexible axis rotation. A:
   Rotating the sketched shape inside the white tubular cavity can
   create artifacts when a fixed rotation axis (dashed) is used. More
   realistically, the shape would revolve like a flexible pipe-cleaner
   (dotted) inside the bearing (gray). B: Fixed rotation
   around an axis :math:`{\mbox{\boldmath ${v}$}}` with a pivot point
   specified by the vector :math:`{\mbox{\boldmath ${u}$}}`.
   C: Subdividing the rotating fragment into slabs with
   separate rotation axes (:math:`\uparrow`) and pivot points
   (:math:`\bullet`) for each slab allows for flexibility. The distance
   between two slabs with indices :math:`n` and :math:`n+1` is
   :math:`\Delta x`.

.. _fig-equipotential:

.. figure:: plots/equipotential.pdf

   Selection of different rotation potentials and definition of
   notation. All four potentials :math:`V` (color coded) are shown for a
   single atom at position :math:`{\mbox{\boldmath ${x}$}}_j(t)`.
   A: Isotropic potential :math:`V^\mathrm{iso}`,
   B: radial motion potential :math:`V^\mathrm{rm}` and
   flexible potential :math:`V^\mathrm{flex}`, C–D: radial
   motion2 potential :math:`V^\mathrm{rm2}` and flexible2 potential
   :math:`V^\mathrm{flex2}` for :math:`\epsilon' = 0`\ nm\ :math:`^2`
   (C) and :math:`\epsilon' = 0.01`\ nm\ :math:`^2`
   (D). The rotation axis is perpendicular to the plane
   and marked by :math:`\otimes`. The light gray contours indicate
   Boltzmann factors :math:`e^{-V/(k_B T)}` in the
   :math:`{\mbox{\boldmath ${x}$}}_j`-plane for :math:`T=300`\ K and
   :math:`k=200`\ kJ/(mol\ :math:`\cdot`\ nm\ :math:`^2`). The green
   arrow shows the direction of the force
   :math:`{\mbox{\boldmath ${F}$}}_{\!j}` acting on atom :math:`j`; the
   blue dashed line indicates the motion of the reference position.

Fixed Axis Rotation
^^^^^^^^^^^^^^^^^^^

Stationary Axis with an Isotropic Potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the fixed axis approach (see :numref:`Fig. (%s) B <fig-rotation>`),
torque on a group of :math:`N` atoms with positions
:math:`{\mbox{\boldmath ${x}$}}_i` (denoted “rotation group”) is applied
by rotating a reference set of atomic positions – usually their initial
positions :math:`{\mbox{\boldmath ${y}$}}_i^0` – at a constant angular
velocity :math:`\omega` around an axis defined by a direction vector
:math:`\hat{{\mbox{\boldmath ${v}$}}}` and a pivot point
:math:`{\mbox{\boldmath ${u}$}}`. To that aim, each atom with
position :math:`{\mbox{\boldmath ${x}$}}_i` is attracted by a “virtual
spring” potential to its moving reference position
:math:`{\mbox{\boldmath ${y}$}}_i = \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})`,
where :math:`\mathbf{\Omega}(t)` is a matrix that describes the rotation
around the axis. In the simplest case, the “springs” are described by a
harmonic potential,

.. math:: V^\mathrm{iso} = \frac{k}{2} \sum_{i=1}^{N} w_i \left[ \mathbf{\Omega}(t)
          ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})  \right]^2
          :label: eqnpotiso

with optional mass-weighted prefactors :math:`w_i = N \, m_i/M` with
total mass :math:`M = \sum_{i=1}^N m_i`. The rotation matrix
:math:`\mathbf{\Omega}(t)` is

.. math::

   \mathbf{\Omega}(t) =  
   \left(   
   \begin{array}{ccc}
   \cos\omega t + v_x^2{\,\xi\,}& v_x v_y{\,\xi\,}- v_z\sin\omega t  & v_x v_z{\,\xi\,}+ v_y\sin\omega t\\
   v_x v_y{\,\xi\,}+ v_z\sin\omega t  & \cos\omega t + v_y^2{\,\xi\,}& v_y v_z{\,\xi\,}- v_x\sin\omega t\\
   v_x v_z{\,\xi\,}- v_y\sin\omega t  & v_y v_z{\,\xi\,}+ v_x\sin\omega t  & \cos\omega t + v_z^2{\,\xi\,}\\
   \end{array}
   \right)

where :math:`v_x`, :math:`v_y`, and :math:`v_z` are the components of
the normalized rotation vector :math:`\hat{{\mbox{\boldmath ${v}$}}}`,
and :math:`{\,\xi\,}:= 1-\cos(\omega t)`. As illustrated in
:numref:`Fig.  (%s) A <fig-equipotential>` for a single atom :math:`j`,
the rotation matrix :math:`\mathbf{\Omega}(t)` operates on the initial
reference positions
:math:`{\mbox{\boldmath ${y}$}}_j^0 = {\mbox{\boldmath ${x}$}}_j(t_0)`
of atom :math:`j` at :math:`t=t_0`. At a later time :math:`t`, the
reference position has rotated away from its initial place (along the
blue dashed line), resulting in the force

.. math:: {\mbox{\boldmath ${F}$}}_{\!j}^\mathrm{iso} 
          = -\nabla_{\!j} \, V^\mathrm{iso} 
          = k \, w_j \left[
          \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_j^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${u}$}} ) \right]
          :label: eqnforcefixed

which is directed towards the reference position.

Pivot-Free Isotropic Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of a fixed pivot vector :math:`{\mbox{\boldmath ${u}$}}` this
potential uses the center of mass :math:`{\mbox{\boldmath ${x}$}}_c` of
the rotation group as pivot for the rotation axis,

.. math:: {\mbox{\boldmath ${x}$}}_c   = \frac{1}{M} \sum_{i=1}^N m_i {\mbox{\boldmath ${x}$}}_i 
          \mbox{\hspace{4ex}and\hspace{4ex}}
          {\mbox{\boldmath ${y}$}}_c^0 = \frac{1}{M} \sum_{i=1}^N m_i {\mbox{\boldmath ${y}$}}_i^0 \ ,
          :label: eqncom

which yields the “pivot-free” isotropic potential

.. math:: \mathchardef\mhyphen="2D
          V^\mathrm{iso\mhyphen pf} = \frac{k}{2} \sum_{i=1}^{N} w_i \left[ \mathbf{\Omega}(t)
          ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c) \right]^2 ,
          :label: eqnpotisopf

with forces

.. math:: \mathchardef\mhyphen="2D
          \mathbf{F}_{\!j}^\mathrm{iso\mhyphen pf} = k \, w_j 
          \left[ 
          \mathbf{\Omega}(t) ( {\mbox{\boldmath ${y}$}}_j^0 - {\mbox{\boldmath ${y}$}}_c^0) 
                           - ( {\mbox{\boldmath ${x}$}}_j   - {\mbox{\boldmath ${x}$}}_c )
          \right] .
          :label: eqnforceisopf

Without mass-weighting, the pivot :math:`{\mbox{\boldmath ${x}$}}_c` is
the geometrical center of the group. 

Parallel Motion Potential Variant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The forces generated by the isotropic potentials
(:eq:`eqns. (%s) <eqnpotiso>` and
:eq:`(%s) <eqnpotisopf>`) also contain components parallel to the
rotation axis and thereby restrain motions along the axis of either the
whole rotation group (in case of :math:`V^\mathrm{iso}`) or within the
rotation group, in case of 

.. math:: 
        \mathchardef\mhyphen="2D
        V^\mathrm{iso\mhyphen pf}
        
For cases where
unrestrained motion along the axis is preferred, we have implemented a
“parallel motion” variant by eliminating all components parallel to the
rotation axis for the potential. This is achieved by projecting the
distance vectors between reference and actual positions

.. math:: {\mbox{\boldmath ${r}$}}_i = \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})

onto the plane perpendicular to the rotation vector,

.. math:: {\mbox{\boldmath ${r}$}}_i^\perp :=  {\mbox{\boldmath ${r}$}}_i - ({\mbox{\boldmath ${r}$}}_i \cdot \hat{{\mbox{\boldmath ${v}$}}})\hat{{\mbox{\boldmath ${v}$}}}
          :label: eqnproject

yielding

.. math:: \begin{aligned}
          \nonumber
          V^\mathrm{pm} &=& \frac{k}{2} \sum_{i=1}^{N} w_i ( {\mbox{\boldmath ${r}$}}_i^\perp )^2 \\
                  &=& \frac{k}{2} \sum_{i=1}^{N} w_i
           \left\lbrace
           \mathbf{\Omega}(t)
             ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})  \right. \nonumber \\
          && \left. - \left\lbrace
          \left[ \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}) \right] \cdot\hat{{\mbox{\boldmath ${v}$}}}
            \right\rbrace\hat{{\mbox{\boldmath ${v}$}}} \right\rbrace^2
          \end{aligned}
          :label: eqnpotpm

and similarly

.. math:: {\mbox{\boldmath ${F}$}}_{\!j}^\mathrm{pm} = k \, w_j \, {\mbox{\boldmath ${r}$}}_j^\perp
          :label: eqnforcepm

Pivot-Free Parallel Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Replacing in :eq:`eqn. (%s) <eqnpotpm>` the fixed pivot
:math:`{\mbox{\boldmath ${u}$}}` by the center of mass
:math:`{\mbox{\boldmath ${x_c}$}}` yields the pivot-free variant of the
parallel motion potential. With

.. math:: 

    \mathchardef\mhyphen="2D
    {\mbox{\boldmath ${s}$}}_i = \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c)

the respective potential and forces are

.. math:: \begin{aligned}
          \mathchardef\mhyphen="2D
          V^\mathrm{pm\mhyphen pf} &=& \frac{k}{2} \sum_{i=1}^{N} w_i ( {\mbox{\boldmath ${s}$}}_i^\perp )^2 \end{aligned}
          :label: eqnpotpmpf

.. math:: \begin{aligned}       
          \mathchardef\mhyphen="2D
          {\mbox{\boldmath ${F}$}}_{\!j}^\mathrm{pm\mhyphen pf} &=& k \, w_j \, {\mbox{\boldmath ${s}$}}_j^\perp
          \end{aligned}
          :label: eqnforcepmpf

Radial Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^

In the above variants, the minimum of the rotation potential is either a
single point at the reference position
:math:`{\mbox{\boldmath ${y}$}}_i` (for the isotropic potentials) or a
single line through :math:`{\mbox{\boldmath ${y}$}}_i` parallel to the
rotation axis (for the parallel motion potentials). As a result, radial
forces restrict radial motions of the atoms. The two subsequent types of
rotation potentials, :math:`V^\mathrm{rm}` and :math:`V^\mathrm{rm2}`, drastically
reduce or even eliminate this effect. The first variant, :math:`V^\mathrm{rm}`
(:numref:`Fig. (%s) B <fig-equipotential>`), eliminates all force
components parallel to the vector connecting the reference atom and the
rotation axis,

.. math:: V^\mathrm{rm} = \frac{k}{2} \sum_{i=1}^{N} w_i \left[
          {\mbox{\boldmath ${p}$}}_i
          \cdot({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}) \right]^2 ,
          :label: eqnpotrm

with

.. math::

   {\mbox{\boldmath ${p}$}}_i := 
   \frac{\hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})} {\| \hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})\|} \ .

This variant depends only on the distance
:math:`{\mbox{\boldmath ${p}$}}_i \cdot ({\mbox{\boldmath ${x}$}}_i -
{\mbox{\boldmath ${u}$}})` of atom :math:`i` from the plane spanned by
:math:`\hat{{\mbox{\boldmath ${v}$}}}` and
:math:`\mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})`.
The resulting force is

.. math:: \mathbf{F}_{\!j}^\mathrm{rm} =
           -k \, w_j \left[ {\mbox{\boldmath ${p}$}}_j\cdot({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${u}$}}) \right] \,{\mbox{\boldmath ${p}$}}_j \,  .
          :label: eqnpotrmforce

Pivot-Free Radial Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Proceeding similar to the pivot-free isotropic potential yields a
pivot-free version of the above potential. With

.. math::

   {\mbox{\boldmath ${q}$}}_i := 
   \frac{\hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0)} {\| \hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0)\|} \, ,

the potential and force for the pivot-free variant of the radial motion
potential read

.. math:: \begin{aligned}
          \mathchardef\mhyphen="2D
          V^\mathrm{rm\mhyphen pf} & = & \frac{k}{2} \sum_{i=1}^{N} w_i \left[
          {\mbox{\boldmath ${q}$}}_i
          \cdot({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c)
          \right]^2 \, , \end{aligned}
          :label: eqnpotrmpf

.. math:: \begin{aligned}       
          \mathchardef\mhyphen="2D
          \mathbf{F}_{\!j}^\mathrm{rm\mhyphen pf} & = &
           -k \, w_j \left[ {\mbox{\boldmath ${q}$}}_j\cdot({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${x}$}}_c) \right] \,{\mbox{\boldmath ${q}$}}_j 
           + k   \frac{m_j}{M} \sum_{i=1}^{N} w_i \left[
           {\mbox{\boldmath ${q}$}}_i\cdot({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c) \right]\,{\mbox{\boldmath ${q}$}}_i \, .
          \end{aligned}
          :label: eqnpotrmpfforce

Radial Motion 2 Alternative Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As seen in :numref:`Fig. (%s) B <fig-equipotential>`, the force
resulting from :math:`V^\mathrm{rm}` still contains a small, second-order
radial component. In most cases, this perturbation is tolerable; if not,
the following alternative, :math:`V^\mathrm{rm2}`, fully eliminates the
radial contribution to the force, as depicted in
:numref:`Fig. (%s) C <fig-equipotential>`,

.. _math-629:

.. math:: V^\mathrm{rm2} = 
          \frac{k}{2} \sum_{i=1}^{N} w_i\, 
          \frac{\left[ (\hat{{\mbox{\boldmath ${v}$}}} \times ( {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}} ))
          \cdot \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) \right]^2}
          {\| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}) \|^2 +
          \epsilon'} \, ,
          :label: eqnpotrm2

where a small parameter :math:`\epsilon'` has been introduced to avoid
singularities. For :math:`\epsilon'=0`\ nm\ :math:`^2`, the
equipotential planes are spanned by :math:`{\mbox{\boldmath ${x}$}}_i -
{\mbox{\boldmath ${u}$}}` and :math:`\hat{{\mbox{\boldmath ${v}$}}}`,
yielding a force perpendicular to
:math:`{\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}`, thus not
contracting or expanding structural parts that moved away from or toward
the rotation axis.

Choosing a small positive :math:`\epsilon'` (*e.g.*,
:math:`\epsilon'=0.01`\ nm\ :math:`^2`,
:numref:`Fig. (%s) D <fig-equipotential>`) in the denominator of
:eq:`eqn. (%s) <eqnpotrm2>` yields a well-defined potential and
continuous forces also close to the rotation axis, which is not the case
for :math:`\epsilon'=0`\ nm\ :math:`^2`
(:numref:`Fig. (%s) C <fig-equipotential>`). With

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${r}$}}_i & := & \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})\\
   {\mbox{\boldmath ${s}$}}_i & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i -
   {\mbox{\boldmath ${u}$}} ) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})
   \| } \equiv \; \Psi_{i} \;\; {\hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}} ) }\\
   \Psi_i^{*}   & := & \frac{1}{ \| \hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}}) \|^2 + \epsilon'}\end{aligned}

the force on atom :math:`j` reads

.. math:: {\mbox{\boldmath ${F}$}}_{\!j}^\mathrm{rm2}  = 
          - k\; 
          \left\lbrace w_j\;
          ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j})\;
          \left[ \frac{\Psi_{\!j}^*   }{\Psi_{\!j}  }  {\mbox{\boldmath ${r}$}}_{\!j} 
               - \frac{\Psi_{\!j}^{ * 2}}{\Psi_{\!j}^3}
               ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j}){\mbox{\boldmath ${s}$}}_j \right]
          \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}} .
          :label: eqnpotrm2force

Pivot-Free Radial Motion 2 Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pivot-free variant of the above potential is

.. math:: \mathchardef\mhyphen="2D
          V{^\mathrm{rm2\mhyphen pf}}= 
          \frac{k}{2} \sum_{i=1}^{N} w_i\, 
          \frac{\left[ (\hat{{\mbox{\boldmath ${v}$}}} \times ( {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c ))
          \cdot \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c) \right]^2}
          {\| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c) \|^2 +
          \epsilon'} \, .
          :label: eqnpotrm2pf

With

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${r}$}}_i & := & \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c)\\
   {\mbox{\boldmath ${s}$}}_i & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i -
   {\mbox{\boldmath ${x}$}}_c ) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c)
   \| } \equiv \; \Psi_{i} \;\; {\hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c ) }\\ \Psi_i^{*}   & := & \frac{1}{ \| \hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c) \|^2 + \epsilon'}\end{aligned}

the force on atom :math:`j` reads

.. math:: \begin{aligned}
          \nonumber
          \mathchardef\mhyphen="2D
          {\mbox{\boldmath ${F}$}}_{\!j}{^\mathrm{rm2\mhyphen pf}}& = &
          - k\; 
          \left\lbrace w_j\;
          ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j})\;
          \left[ \frac{\Psi_{\!j}^*   }{\Psi_{\!j}  } {\mbox{\boldmath ${r}$}}_{\!j} 
               - \frac{\Psi_{\!j}^{ * 2}}{\Psi_{\!j}^3}
               ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j}){\mbox{\boldmath ${s}$}}_j \right]
          \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}}\\
               & &
          + k\;\frac{m_j}{M} \left\lbrace \sum_{i=1}^{N}
          w_i\;({\mbox{\boldmath ${s}$}}_i\cdot{\mbox{\boldmath ${r}$}}_i) \; 
          \left[ \frac{\Psi_i^*   }{\Psi_i  }  {\mbox{\boldmath ${r}$}}_i
               - \frac{\Psi_i^{ * 2}}{\Psi_i^3} ({\mbox{\boldmath ${s}$}}_i\cdot{\mbox{\boldmath ${r}$}}_i )\;
               {\mbox{\boldmath ${s}$}}_i \right] \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}} \, .
          \end{aligned}
          :label: eqnpotrm2pfforce

Flexible Axis Rotation
~~~~~~~~~~~~~~~~~~~~~~

As sketched in :numref:`Fig. (%s) <fig-rotation>` A–B, the rigid body
behavior of the fixed axis rotation scheme is a drawback for many
applications. In particular, deformations of the rotation group are
suppressed when the equilibrium atom positions directly depend on the
reference positions. To avoid this limitation,
:eq:`eqns. (%s) <eqnpotrmpf>` and :eq:`(%s) <eqnpotrm2pf>`
will now be generalized towards a “flexible axis” as sketched in
:numref:`Fig. (%s) C <fig-rotation>`. This will be achieved by
subdividing the rotation group into a set of equidistant slabs
perpendicular to the rotation vector, and by applying a separate
rotation potential to each of these slabs.
:numref:`Fig. (%s) C <fig-rotation>` shows the midplanes of the slabs
as dotted straight lines and the centers as thick black dots.

To avoid discontinuities in the potential and in the forces, we define
“soft slabs” by weighing the contributions of each slab :math:`n` to the
total potential function :math:`V^\mathrm{flex}` by a Gaussian function

.. _math-639:

.. math:: g_n({\mbox{\boldmath ${x}$}}_i) = \Gamma \ \mbox{exp} \left(
          -\frac{\beta_n^2({\mbox{\boldmath ${x}$}}_i)}{2\sigma^2}  \right) ,
          :label: eqngaussian

centered at the midplane of the :math:`n`\ th slab. Here :math:`\sigma`
is the width of the Gaussian function, :math:`\Delta x` the distance
between adjacent slabs, and

.. math:: \beta_n({\mbox{\boldmath ${x}$}}_i) := {\mbox{\boldmath ${x}$}}_i \cdot \hat{{\mbox{\boldmath ${v}$}}} - n \, \Delta x \, .

.. _fig-gaussian:

.. figure:: plots/gaussians.pdf

   Gaussian functions :math:`g_n` centered at :math:`n \, \Delta x` for
   a slab distance :math:`\Delta x = 1.5` nm and :math:`n \geq -2`.
   Gaussian function :math:`g_0` is highlighted in bold; the dashed line
   depicts the sum of the shown Gaussian functions.

A most convenient choice is :math:`\sigma = 0.7 \Delta x` and

.. math::

   1/\Gamma = \sum_{n \in Z}
   \mbox{exp}
   \left(-\frac{(n - \frac{1}{4})^2}{2\cdot 0.7^2}\right)
   \approx 1.75464 \, ,

which yields a nearly constant sum, essentially independent of
:math:`{\mbox{\boldmath ${x}$}}_i` (dashed line in
:numref:`Fig. (%s) <fig-gaussian>`), *i.e.*,

.. math:: \sum_{n \in Z} g_n({\mbox{\boldmath ${x}$}}_i) =  1 + \epsilon({\mbox{\boldmath ${x}$}}_i) \, ,
          :label: eqnnormal

with
:math:`| \epsilon({\mbox{\boldmath ${x}$}}_i) | < 1.3\cdot 10^{-4}`.
This choice also implies that the individual contributions to the force
from the slabs add up to unity such that no further normalization is
required.

To each slab center :math:`{\mbox{\boldmath ${x}$}}_c^n`, all atoms
contribute by their Gaussian-weighted (optionally also mass-weighted)
position vectors
:math:`g_n({\mbox{\boldmath ${x}$}}_i) \, {\mbox{\boldmath ${x}$}}_i`.
The instantaneous slab centers :math:`{\mbox{\boldmath ${x}$}}_c^n` are
calculated from the current positions
:math:`{\mbox{\boldmath ${x}$}}_i`,

.. math::  {\mbox{\boldmath ${x}$}}_c^n =
           \frac{\sum_{i=1}^N g_n({\mbox{\boldmath ${x}$}}_i) \, m_i \, {\mbox{\boldmath ${x}$}}_i}
                {\sum_{i=1}^N g_n({\mbox{\boldmath ${x}$}}_i) \, m_i} \, ,\\
           :label: eqndefx0 

while the reference centers :math:`{\mbox{\boldmath ${y}$}}_c^n` are
calculated from the reference positions
:math:`{\mbox{\boldmath ${y}$}}_i^0`,

.. math:: {\mbox{\boldmath ${y}$}}_c^n =
          \frac{\sum_{i=1}^N g_n({\mbox{\boldmath ${y}$}}_i^0) \, m_i \, {\mbox{\boldmath ${y}$}}_i^0}
               {\sum_{i=1}^N g_n({\mbox{\boldmath ${y}$}}_i^0) \, m_i} \, .
          :label: eqndefy0

Due to the rapid decay of :math:`g_n`, each slab will essentially
involve contributions from atoms located within :math:`\approx
3\Delta x` from the slab center only.

Flexible Axis Potential
^^^^^^^^^^^^^^^^^^^^^^^

We consider two flexible axis variants. For the first variant, the slab
segmentation procedure with Gaussian weighting is applied to the radial
motion potential
(:eq:`eqn. (%s) <eqnpotrmpf>` / :numref:`Fig. (%s) B <fig-equipotential>`),
yielding as the contribution of slab :math:`n`

.. math::  V^n = 
           \frac{k}{2} \sum_{i=1}^{N} w_i \, g_n({\mbox{\boldmath ${x}$}}_i) 
           \left[
           {\mbox{\boldmath ${q}$}}_i^n
           \cdot
            ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n) 
           \right]^2  ,
           :label: eqnflexpot

and a total potential function

.. math:: V^\mathrm{flex} = \sum_n V^n \, .
          :label: eqnpotflex

Note that the global center of mass :math:`{\mbox{\boldmath ${x}$}}_c`
used in :eq:`eqn. (%s) <eqnpotrmpf>` is now replaced by
:math:`{\mbox{\boldmath ${x}$}}_c^n`, the center of mass of the slab.
With

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${q}$}}_i^n & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times
   \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n) }{ \| \hat{{\mbox{\boldmath ${v}$}}}
   \times \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n) \| } \\
   b_i^n         & := & {\mbox{\boldmath ${q}$}}_i^n \cdot ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n) \, ,\end{aligned}

the resulting force on atom :math:`j` reads

.. math:: \begin{aligned}
          \nonumber\hspace{-15mm}
          {\mbox{\boldmath ${F}$}}_{\!j}^\mathrm{flex} &=&
          - \, k \, w_j \sum_n g_n({\mbox{\boldmath ${x}$}}_j) \, b_j^n \left\lbrace  {\mbox{\boldmath ${q}$}}_j^n -
          b_j^n \frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{2\sigma^2} \hat{{\mbox{\boldmath ${v}$}}} \right\rbrace \\ & &
          + \, k \, m_j \sum_n \frac{g_n({\mbox{\boldmath ${x}$}}_j)}{\sum_h g_n({\mbox{\boldmath ${x}$}}_h)}
          \sum_{i=1}^{N} w_i \, g_n({\mbox{\boldmath ${x}$}}_i) \, b_i^n \left\lbrace 
          {\mbox{\boldmath ${q}$}}_i^n -\frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{\sigma^2}
          \left[ {\mbox{\boldmath ${q}$}}_i^n \cdot ({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${x}$}}_c^n )\right]
          \hat{{\mbox{\boldmath ${v}$}}} \right\rbrace .
          \end{aligned}
          :label: eqnpotflexforce

Note that for :math:`V^\mathrm{flex}`, as defined, the slabs are fixed in
space and so are the reference centers
:math:`{\mbox{\boldmath ${y}$}}_c^n`. If during the simulation the
rotation group moves too far in :math:`{\mbox{\boldmath ${v}$}}`
direction, it may enter a region where – due to the lack of nearby
reference positions – no reference slab centers are defined, rendering
the potential evaluation impossible. We therefore have included a
slightly modified version of this potential that avoids this problem by
attaching the midplane of slab :math:`n=0` to the center of mass of the
rotation group, yielding slabs that move with the rotation group. This
is achieved by subtracting the center of mass
:math:`{\mbox{\boldmath ${x}$}}_c` of the group from the positions,

.. math:: \tilde{{\mbox{\boldmath ${x}$}}}_i = {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c \, , \mbox{\ \ \ and \ \ } 
          \tilde{{\mbox{\boldmath ${y}$}}}_i^0 = {\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0 \, ,
          :label: eqntrafo

such that

.. math:: \begin{aligned}
          \mathchardef\mhyphen="2D
          V^\mathrm{flex\mhyphen t} 
            & = & \frac{k}{2} \sum_n \sum_{i=1}^{N} w_i \, g_n(\tilde{{\mbox{\boldmath ${x}$}}}_i)
            \left[ \frac{\hat{{\mbox{\boldmath ${v}$}}} \times \mathbf{\Omega}(t)(\tilde{{\mbox{\boldmath ${y}$}}}_i^0
            - \tilde{{\mbox{\boldmath ${y}$}}}_c^n) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times
          \mathbf{\Omega}(t)(\tilde{{\mbox{\boldmath ${y}$}}}_i^0 -
          \tilde{{\mbox{\boldmath ${y}$}}}_c^n) \| }
          \cdot
           (\tilde{{\mbox{\boldmath ${x}$}}}_i - \tilde{{\mbox{\boldmath ${x}$}}}_c^n) 
          \right]^2 .
          \end{aligned}
          :label: eqnpotflext

To simplify the force derivation, and for efficiency reasons, we here
assume :math:`{\mbox{\boldmath ${x}$}}_c` to be constant, and thus
:math:`\partial {\mbox{\boldmath ${x}$}}_c / \partial x =
\partial {\mbox{\boldmath ${x}$}}_c / \partial y = \partial {\mbox{\boldmath ${x}$}}_c / \partial z = 0`.
The resulting force error is small (of order :math:`O(1/N)` or
:math:`O(m_j/M)` if mass-weighting is applied) and can therefore be
tolerated. With this assumption, the forces

.. math::
    \mathchardef\mhyphen="2D  
    {\mbox{\boldmath ${F}$}}^\mathrm{flex\mhyphen t}
   
have the same form as
:eq:`eqn. (%s) <eqnpotflexforce>`.

Flexible Axis 2 Alternative Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this second variant, slab segmentation is applied to
:math:`V^\mathrm{rm2}` (:eq:`eqn. (%s) <eqnpotrm2pf>`), resulting in
a flexible axis potential without radial force contributions
(:numref:`Fig. (%s) C <fig-equipotential>`),

.. math::   V{^\mathrm{flex2}}= 
            \frac{k}{2} \sum_{i=1}^{N} \sum_n w_i\,g_n({\mbox{\boldmath ${x}$}}_i) 
            \frac{\left[ (\hat{{\mbox{\boldmath ${v}$}}} \times ( {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n ))
            \cdot \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n) \right]^2}
            {\| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n) \|^2 +
            \epsilon'} \, .
            :label: eqnpotflex2

With

.. math:: \begin{aligned}
          {\mbox{\boldmath ${r}$}}_i^n & := & \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n)\\
          {\mbox{\boldmath ${s}$}}_i^n & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i -
          {\mbox{\boldmath ${x}$}}_c^n ) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n)
          \| } \equiv \; \psi_{i} \;\; {\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c^n ) }\\
          \psi_i^{*}     & := & \frac{1}{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c^n) \|^2 + \epsilon'}\\
          W_j^n          & := & \frac{g_n({\mbox{\boldmath ${x}$}}_j)\,m_j}{\sum_h g_n({\mbox{\boldmath ${x}$}}_h)\,m_h}\\
          {\mbox{\boldmath ${S}$}}^n   & := & 
          \sum_{i=1}^{N} w_i\;g_n({\mbox{\boldmath ${x}$}}_i)
          \; ({\mbox{\boldmath ${s}$}}_i^n\cdot{\mbox{\boldmath ${r}$}}_i^n)
          \left[ \frac{\psi_i^*   }{\psi_i  }  {\mbox{\boldmath ${r}$}}_i^n
               - \frac{\psi_i^{ * 2}}{\psi_i^3} ({\mbox{\boldmath ${s}$}}_i^n\cdot{\mbox{\boldmath ${r}$}}_i^n )\;
               {\mbox{\boldmath ${s}$}}_i^n \right] 
          \end{aligned}
          :label: eqnSn

the force on atom :math:`j` reads

.. math:: \begin{aligned}
          \nonumber
          {\mbox{\boldmath ${F}$}}_{\!j}{^\mathrm{flex2}}& = &
          - k\; 
          \left\lbrace \sum_n w_j\;g_n({\mbox{\boldmath ${x}$}}_j)\;
          ({\mbox{\boldmath ${s}$}}_j^n\cdot{\mbox{\boldmath ${r}$}}_{\!j}^n)\;
          \left[ \frac{\psi_j^*   }{\psi_j  }  {\mbox{\boldmath ${r}$}}_{\!j}^n 
               - \frac{\psi_j^{ * 2}}{\psi_j^3} ({\mbox{\boldmath ${s}$}}_j^n\cdot{\mbox{\boldmath ${r}$}}_{\!j}^n)\;
               {\mbox{\boldmath ${s}$}}_{\!j}^n \right] \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}} \\
          \nonumber
          & &
          + k \left\lbrace \sum_n W_{\!j}^n \, {\mbox{\boldmath ${S}$}}^n \right\rbrace \times
          \hat{{\mbox{\boldmath ${v}$}}}
          - k \left\lbrace \sum_n W_{\!j}^n \; \frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{\sigma^2} \frac{1}{\psi_j}\;\; 
          {\mbox{\boldmath ${s}$}}_j^n \cdot 
          {\mbox{\boldmath ${S}$}}^n \right\rbrace \hat{{\mbox{\boldmath ${v}$}}}\\ 
          & & 
          + \frac{k}{2} \left\lbrace \sum_n w_j\;g_n({\mbox{\boldmath ${x}$}}_j)
          \frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{\sigma^2} 
          \frac{\psi_j^*}{\psi_j^2}( {\mbox{\boldmath ${s}$}}_j^n \cdot {\mbox{\boldmath ${r}$}}_{\!j}^n )^2 \right\rbrace
          \hat{{\mbox{\boldmath ${v}$}}} .
          \end{aligned}
          :label: eqnpotflex2force

Applying transformation (:eq:`(%s) <eqntrafo>`) yields a
“translation-tolerant” version of the flexible2 potential,

.. math::
    \mathchardef\mhyphen="2D
    V{^\mathrm{flex2\mhyphen t}}
    
Again, assuming that
:math:`\partial {\mbox{\boldmath ${x}$}}_c / \partial x`,
:math:`\partial {\mbox{\boldmath ${x}$}}_c /
\partial y`, :math:`\partial {\mbox{\boldmath ${x}$}}_c / \partial z`
are small, the resulting equations for

.. math::
    \mathchardef\mhyphen="2D
    V{^\mathrm{flex2\mhyphen t}}
    
and

.. math::
    \mathchardef\mhyphen="2D
    {\mbox{\boldmath ${F}$}}{^\mathrm{flex2\mhyphen t}}
   
are similar
to those of :math:`V^\mathrm{flex2}` and
:math:`{\mbox{\boldmath ${F}$}}^\mathrm{flex2}`.

Usage
~~~~~

To apply enforced rotation, the particles :math:`i` that are to be
subjected to one of the rotation potentials are defined via index groups
``rot-group0``, ``rot-group1``, etc., in the
:ref:`mdp` input file. The reference positions
:math:`{\mbox{\boldmath ${y}$}}_i^0` are read from a special
:ref:`trr` file provided to :ref:`grompp <gmx grompp>`. If no such
file is found, :math:`{\mbox{\boldmath ${x}$}}_i(t=0)` are used as
reference positions and written to :ref:`trr` such that they
can be used for subsequent setups. All parameters of the potentials such
as :math:`k`, :math:`\epsilon'`, etc.
(:numref:`Table (%s) <tab-vars>`) are provided as :ref:`mdp`
parameters; ``rot-type`` selects the type of the potential.
The option ``rot-massw`` allows to choose whether or not to
use mass-weighted averaging. For the flexible potentials, a cutoff value
:math:`g_n^\mathrm{min}` (typically :math:`g_n^\mathrm{min}=0.001`)
makes shure that only significant contributions to :math:`V` and
:math:`{\mbox{\boldmath ${F}$}}` are evaluated, *i.e.* terms with
:math:`g_n({\mbox{\boldmath ${x}$}}) < g_n^\mathrm{min}` are omitted.
:numref:`Table (%s) <tab-quantities>` summarizes observables that are
written to additional output files and which are described below.

.. |ROTISO| replace:: V\ :math:`^{\mathrm{iso}}`
.. |ROTISOPF| replace:: V\ :math:`^{\mathrm{iso--pf}}`
.. |ROTPM| replace:: V\ :math:`^{\mathrm{pm}}`
.. |ROTPMPF| replace:: V\ :math:`^{\mathrm{pm--pf}}`
.. |ROTRM| replace:: V\ :math:`^{\mathrm{rm}}`
.. |ROTRMPF| replace:: V\ :math:`^{\mathrm{rm--pf}}`
.. |ROTRM2| replace:: V\ :math:`^{\mathrm{rm2}}`
.. |ROTRM2PF| replace:: V\ :math:`^{\mathrm{rm2--pf}}`
.. |ROTFL| replace:: V\ :math:`^{\mathrm{flex}}`
.. |ROTFLT| replace:: V\ :math:`^{\mathrm{flex--t}}`
.. |ROTFL2| replace:: V\ :math:`^{\mathrm{flex2}}`
.. |ROTFLT2| replace:: V\ :math:`^{\mathrm{flex2--t}}`
.. |KUNIT| replace:: :math:`\frac{\mathrm{kJ}}{\mathrm{mol} \cdot \mathrm{nm}^2}`
.. |BFX| replace:: **x**
.. |KMA| replace:: :math:`k`
.. |VECV| replace:: :math:`\hat{\mbox{boldmath{$v$}}}`
.. |VECU| replace:: :math:`{\mbox{\boldmath{$u$}}}`
.. |OMEG| replace:: :math:`\omega`
.. |EPSP| replace:: :math:`{\epsilon}'`
.. |DELX| replace:: :math:`{\Delta}x`
.. |GMIN| replace:: :math:`g_n^\mathrm{min}`
.. |CIPS| replace:: :math:`^\circ`\ /ps
.. |NM2| replace:: nm\ :math:`^2`
.. |REF1| replace:: (\ :eq:`eqnpotiso`\ )
.. |REF2| replace:: (\ :eq:`eqnpotisopf`\ )
.. |REF3| replace:: (\ :eq:`eqnpotpm`\ )
.. |REF4| replace:: (\ :eq:`eqnpotpmpf`\ )
.. |REF5| replace:: (\ :eq:`eqnpotrm`\ )
.. |REF6| replace:: (\ :eq:`eqnpotrmpf`\ )
.. |REF7| replace:: (\ :eq:`eqnpotrm2`\ )
.. |REF8| replace:: (\ :eq:`eqnpotrm2pf`\ )
.. |REF9| replace:: (\ :eq:`eqnpotflex`\ )
.. |REF10| replace:: (\ :eq:`eqnpotflext`\ )
.. |REF11| replace:: (\ :eq:`eqnpotflex2`\ )

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


.. |VT|      replace:: :math:`V(t)`
.. |THET|    replace:: :math:`\theta_\mathrm{ref}(t)`
.. |THETAV|  replace:: :math:`\theta_\mathrm{av}(t)`
.. |THETFIT| replace:: :math:`\theta_\mathrm{fit}(t)`, :math:`\theta_\mathrm{fit}(t,n)`
.. |YVEC|    replace:: :math:`{\mbox{boldmath{$y$}}}_{0}(n)`, :math:`{\mbox{boldmath{$x$}}}_{0}(t,n)`
.. |TAUT|    replace:: :math:`\tau(t)`
.. |TAUTN|   replace:: :math:`\tau(t,n)`
.. |REFT|  replace:: :numref:`see Table (%s) <tab-vars>`
.. |REFEQ| replace:: :math:`\theta_\mathrm{ref}(t)=\omega t`
.. |REF12| replace:: (\ :eq:`eqnavangle`\ )               
.. |REF13| replace:: (\ :eq:`eqnrmsdfit`\ )               
.. |REF14| replace:: (\ :eq:`eqndefx0`\ ,\ :eq:`eqndefy0`\ )
.. |REF15| replace:: (\ :eq:`eqntorque`\ )                

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
rotation axis through pivot point :math:`{\mbox{\boldmath ${u}$}}` (see
:eq:`eqn. (%s) <eqnproject>` for the definition of
:math:`\perp`),

.. math::

   \cos \theta_i = 
   \frac{({\mbox{\boldmath ${y}$}}_i-{\mbox{\boldmath ${u}$}})^\perp \cdot ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}})^\perp}
        { \| ({\mbox{\boldmath ${y}$}}_i-{\mbox{\boldmath ${u}$}})^\perp \cdot ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}})^\perp
        \| } \ .

The sign of :math:`\theta_\mathrm{av}` is chosen such that
:math:`\theta_\mathrm{av} > 0` if the actual structure rotates ahead of
the reference.

Angle of Rotation Groups: Flexible Axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For flexible axis rotation, two outputs are provided, the angle of the
entire rotation group, and separate angles for the segments in the
slabs. The angle of the entire rotation group is determined by an RMSD
fit of :math:`{\mbox{\boldmath ${x}$}}_i` to the reference positions
:math:`{\mbox{\boldmath ${y}$}}_i^0` at :math:`t=0`, yielding
:math:`\theta_\mathrm{fit}` as the angle by which the reference has to
be rotated around :math:`\hat{{\mbox{\boldmath ${v}$}}}` for the optimal
fit,

.. math::  \mathrm{RMSD} \big( {\mbox{\boldmath ${x}$}}_i,\ \mathbf{\Omega}(\theta_\mathrm{fit})
           {\mbox{\boldmath ${y}$}}_i^0 \big) \stackrel{!}{=} \mathrm{min} \, .
           :label: eqnrmsdfit

To determine the local angle for each slab :math:`n`, both reference
and actual positions are weighted with the Gaussian function of slab
:math:`n`, and :math:`\theta_\mathrm{fit}(t,n)` is calculated as in
:eq:`eqn. (%s) <eqnrmsdfit>`) from the Gaussian-weighted
positions.

For all angles, the :ref:`mdp` input option
``rot-fit-method`` controls whether a normal RMSD fit is
performed or whether for the fit each position
:math:`{\mbox{\boldmath ${x}$}}_i` is put at the same distance to the
rotation axis as its reference counterpart
:math:`{\mbox{\boldmath ${y}$}}_i^0`. In the latter case, the RMSD
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

The torque :math:`{\mbox{\boldmath ${\tau}$}}(t)` exerted by the
rotation potential is calculated for fixed axis rotation via

.. math:: {\mbox{\boldmath ${\tau}$}}(t) = \sum_{i=1}^{N} {\mbox{\boldmath ${r}$}}_i(t) \times {\mbox{\boldmath ${f}$}}_{\!i}^\perp(t) ,
          :label: eqntorque

where :math:`{\mbox{\boldmath ${r}$}}_i(t)` is the distance vector from
the rotation axis to :math:`{\mbox{\boldmath ${x}$}}_i(t)` and
:math:`{\mbox{\boldmath ${f}$}}_{\!i}^\perp(t)` is the force component
perpendicular to :math:`{\mbox{\boldmath ${r}$}}_i(t)` and
:math:`\hat{{\mbox{\boldmath ${v}$}}}`. For flexible axis rotation,
torques :math:`{\mbox{\boldmath ${\tau}$}}_{\!n}` are calculated for
each slab using the local rotation axis of the slab and the
Gaussian-weighted positions.

Electric fields
---------------

A pulsed and oscillating electric field can be applied according to:

.. math:: E(t) = E_0 \exp\left[-\frac{(t-t_0)^2}{2\sigma^2}\right]\cos\left[\omega (t-t_0)\right]
          :label: eq-efield

where :math:`E_0` is the field strength, the angular frequency
:math:`{\mbox{$\omega = 2\pi c/\lambda$}}`, :math:`t_0` is the time
at of the peak in the field strength and :math:`\sigma` is the with of
the pulse. Special cases occur when :math:`\sigma` = 0 (non-pulsed
field) and for :math:`\omega` is 0 (static field).

This simulated laser-pulse was applied to
simulations of melting ice Caleman and Spoel (2008). A pulsed electric
field may look like :numref:`Fig. (%s) <fig-field>`. In the supporting
information of that paper the impact of an applied electric field on a
system under periodic boundary conditions is analyzed. It is described
that the effective electric field under PBC is larger than the applied
field, by a factor depending on the size of the box and the dielectric
properties of molecules in the box. For a system with static dielectric
properties this factor can be corrected for. But for a system where the
dielectric varies over time, for example a membrane protein with a pore
that opens and closes during the simulation, this way of applying an
electric field is not useful. In such cases one can use the
computational electrophysiology protocol described in the next section
(sec. 
).

.. _fig-field:

.. figure:: plots/field.eps

   A simulated laser pulse in |Gromacs|.

Electric fields are applied when the following options are specified in
the :ref:`grompp <gmx grompp>` :ref:`mdp` file. You specify, in order, :math:`E_0`,
:math:`\omega`, :math:`t_0` and :math:`\sigma`:

::

    ElectricField-x = 0.04 0       0     0

yields a static field with :math:`E_0` = 0.04 V/nm in the X-direction.
In contrast,

::

    ElectricField-x = 2.0  150     5     0

yields an oscillating electric field with :math:`E_0` = 2 V/nm,
:math:`\omega` = 150/ps and :math:`t_0` = 5 ps. Finally

::

    ElectricField-x = 2.0  150     5     1

yields an pulsed-oscillating electric field with :math:`E_0` = 2 V/nm,
:math:`\omega` = 150/ps and :math:`t_0` = 5 ps and :math:`\sigma` = 1
ps. Read more in ref. Caleman and Spoel (2008). Note that the input file
format is changed from the undocumented older version. A figure like
:numref:`Fig. (%s) <fig-field>` may be produced by passing the
``-field`` option to :ref:`gmx mdrun`.

Computational Electrophysiology
-------------------------------

The Computational Electrophysiology (CompEL) protocol Kutzner et al.
(2011) allows the simulation of ion flux through membrane channels,
driven by transmembrane potentials or ion concentration gradients. Just
as in real cells, CompEL establishes transmembrane potentials by
sustaining a small imbalance of charges :math:`\Delta q` across the
membrane, which gives rise to a potential difference :math:`\Delta U`
according to the membrane capacitance:

.. math:: \Delta U = \Delta q / C_{membrane}

The transmembrane electric field and concentration gradients are
controlled by :ref:`mdp` options, which allow the user to set
reference counts for the ions on either side of the membrane. If a
difference between the actual and the reference numbers persists over a
certain time span, specified by the user, a number of ion/water pairs
are exchanged between the compartments until the reference numbers are
restored. Alongside the calculation of channel conductance and ion
selectivity, CompEL simulations also enable determination of the channel
reversal potential, an important characteristic obtained in
electrophysiology experiments.

In a CompEL setup, the simulation system is divided into two
compartments **A** and **B** with independent ion concentrations. This
is best achieved by using double bilayer systems with a copy (or copies)
of the channel/pore of interest in each bilayer
(:numref:`Fig. (%s) <fig-compelsetup>` A, B). If the channel axes
point in the same direction, channel flux is observed simultaneously at
positive and negative potentials in this way, which is for instance
important for studying channel rectification.

.. _fig-compelsetup:

.. figure:: plots/compelsetup.pdf

   Typical double-membrane setup for CompEL simulations (A, B).
   Ion/water molecule exchanges will be performed as needed between the
   two light blue volumes around the dotted black lines (A). Plot (C)
   shows the potential difference :math:`\Delta U` resulting from the
   selected charge imbalance :math:`\Delta q_{ref}` between the
   compartments.

The potential difference :math:`\Delta U` across the membrane is easily
calculated with the :ref:`gmx potential <gmx potential>` utility. By this, the
potential drop along :math:`z` or the pore axis is exactly known in each
time interval of the simulation (Fig. 
C). Type and number of ions :math:`n_i` of charge :math:`q_i`,
traversing the channel in the simulation, are written to the
``swapions.xvg`` output file, from which the average channel
conductance :math:`G` in each interval :math:`\Delta t` is determined
by:

.. math:: G = \frac{\sum_{i} n_{i}q_{i}}{\Delta t \, \Delta U} \, .

The ion selectivity is calculated as the number flux ratio of different
species. Best results are obtained by averaging these values over
several overlapping time intervals.

The calculation of reversal potentials is best achieved using a small
set of simulations in which a given transmembrane concentration gradient
is complemented with small ion imbalances of varying magnitude. For
example, if one compartment contains 1M salt and the other 0.1M, and
given charge neutrality otherwise, a set of simulations with
:math:`\Delta q = 0\,e`, :math:`\Delta q = 2\,e`,
:math:`\Delta q = 4\,e` could be used. Fitting a straight line through
the current-voltage relationship of all obtained :math:`I`-:math:`U`
pairs near zero current will then yield :math:`U_{rev}`.

Usage
^^^^^

The following :ref:`mdp` options control the CompEL protocol:

::

    swapcoords     = Z        ; Swap positions: no, X, Y, Z
    swap-frequency = 100      ; Swap attempt frequency

Choose ``Z`` if your membrane is in the :math:`xy`-plane
(:numref:`Fig. (%s) <fig-compelsetup>`). Ions will be exchanged
between compartments depending on their :math:`z`-positions alone.
``swap-frequency`` determines how often a swap attempt will
be made. This step requires that the positions of the split groups, the
ions, and possibly the solvent molecules are communicated between the
parallel processes, so if chosen too small it can decrease the
simulation performance. The ``Position`` swapping entry in
the cycle and time accounting table at the end of the
``md.log`` file summarizes the amount of runtime spent in
the swap module.

::

    split-group0   = channel0 ; Defines compartment boundary
    split-group1   = channel1 ; Defines other compartment boundary
    massw-split0   = no       ; use mass-weighted center?
    massw-split1   = no

``split-group0`` and ``split-group1`` are two
index groups that define the boundaries between the two compartments,
which are usually the centers of the channels. If
``massw-split0`` or ``massw-split1`` are set to
``yes``, the center of mass of each index group is used as
boundary, here in :math:`z`-direction. Otherwise, the geometrical
centers will be used (:math:`\times` in
:numref:`Fig. (%s) <fig-compelsetup>` A). If, such as here, a membrane
channel is selected as split group, the center of the channel will
define the dividing plane between the compartments (dashed horizontal
lines). All index groups must be defined in the index file.

If, to restore the requested ion counts, an ion from one compartment has
to be exchanged with a water molecule from the other compartment, then
those molecules are swapped which have the largest distance to the
compartment-defining boundaries (dashed horizontal lines). Depending on
the ion concentration, this effectively results in exchanges of
molecules between the light blue volumes. If a channel is very
asymmetric in :math:`z`-direction and would extend into one of the swap
volumes, one can offset the swap exchange plane with the
``bulk-offset`` parameter. A value of 0.0 means no offset
:math:`b`, values :math:`-1.0 < b < 0` move the swap exchange plane
closer to the lower, values :math:`0 < b < 1.0` closer to the upper
membrane. :numref:`Fig. (%s) <fig-compelsetup>` A (left) depicts that
for the **A** compartment.

::

    solvent-group  = SOL      ; Group containing the solvent molecules
    iontypes       = 3        ; Number of different ion types to control
    iontype0-name  = NA       ; Group name of the ion type
    iontype0-in-A  = 51       ; Reference count of ions of type 0 in A
    iontype0-in-B  = 35       ; Reference count of ions of type 0 in B
    iontype1-name  = K
    iontype1-in-A  = 10
    iontype1-in-B  = 38
    iontype2-name  = CL
    iontype2-in-A  = -1
    iontype2-in-B  = -1

The group name of solvent molecules acting as exchange partners for the
ions has to be set with ``solvent-group``. The number of
different ionic species under control of the CompEL protocol is given by
the ``iontypes`` parameter, while
``iontype0-name`` gives the name of the index group
containing the atoms of this ionic species. The reference number of ions
of this type can be set with the ``iontype0-in-A`` and
``iontype0-in-B`` options for compartments **A** and **B**,
respectively. Obviously, the sum of ``iontype0-in-A`` and
``iontype0-in-B`` needs to equal the number of ions in the
group defined by ``iontype0-name``. A reference number of
``-1`` means: use the number of ions as found at the
beginning of the simulation as the reference value.

::

    coupl-steps    = 10       ; Average over these many swap steps
    threshold      = 1        ; Do not swap if < threshold

If ``coupl-steps`` is set to 1, then the momentary ion
distribution determines whether ions are exchanged.
``coupl-steps > 1`` will use the time-average of ion
distributions over the selected number of attempt steps instead. This
can be useful, for example, when ions diffuse near compartment
boundaries, which would lead to numerous unproductive ion exchanges. A
``threshold`` of 1 means that a swap is performed if the
average ion count in a compartment differs by at least 1 from the
requested values. Higher thresholds will lead to toleration of larger
differences. Ions are exchanged until the requested number :math:`\pm`
the threshold is reached.

::

    cyl0-r         = 5.0      ; Split cylinder 0 radius (nm)
    cyl0-up        = 0.75     ; Split cylinder 0 upper extension (nm)
    cyl0-down      = 0.75     ; Split cylinder 0 lower extension (nm)
    cyl1-r         = 5.0      ; same for other channel
    cyl1-up        = 0.75
    cyl1-down      = 0.75

The cylinder options are used to define virtual geometric cylinders
around the channel’s pore to track how many ions of which type have
passed each channel. Ions will be counted as having traveled through a
channel according to the definition of the channel’s cylinder radius,
upper and lower extension, relative to the location of the respective
split group. This will not affect the actual flux or exchange, but will
provide you with the ion permeation numbers across each of the channels.
Note that an ion can only be counted as passing through a particular
channel if it is detected *within* the defined split cylinder in a swap
step. If ``swap-frequency`` is chosen too high, a particular
ion may be detected in compartment **A** in one swap step, and in
compartment **B** in the following swap step, so it will be unclear
through which of the channels it has passed.

A double-layered system for CompEL simulations can be easily prepared by
duplicating an existing membrane/channel MD system in the direction of
the membrane normal (typically :math:`z`) with 
:ref:`gmx editconf` ``-translate 0 0 <l_z>``, where ``l_z`` is the box
length in that direction. If you have already defined index groups for
the channel for the single-layered system, :ref:`gmx make_ndx`
``-n index.ndx -twin`` will provide you with the groups for the
double-layered system.

To suppress large fluctuations of the membranes along the swap
direction, it may be useful to apply a harmonic potential (acting only
in the swap dimension) between each of the two channel and/or bilayer
centers using umbrella pulling (see
section 
).

Multimeric channels
^^^^^^^^^^^^^^^^^^^

If a split group consists of more than one molecule, the correct PBC
image of all molecules with respect to each other has to be chosen such
that the channel center can be correctly determined. |Gromacs| assumes
that the starting structure in the :ref:`tpr` file has the
correct PBC representation. Set the following environment variable to
check whether that is the case:

-  ``GMX_COMPELDUMP``: output the starting structure after
   it has been made whole to :ref:`pdb` file.

Calculating a PMF using the free-energy code
--------------------------------------------

The free-energy
coupling-parameter approach (see sec. 
)
provides several ways to calculate potentials of mean force. A potential
of mean force between two atoms can be calculated by connecting them
with a harmonic potential or a constraint. For this purpose there are
special potentials that avoid the generation of extra exclusions,
see sec. 
. When the position of the minimum
or the constraint length is 1 nm more in state B than in state A, the
restraint or constraint force is given by
:math:`\partial H/\partial \lambda`. The distance between the atoms can
be changed as a function of :math:`\lambda` and time by setting
``delta-lambda`` in the :ref:`mdp` file. The
results should be identical (although not numerically due to the
different implementations) to the results of the pull code with umbrella
sampling and constraint pulling. Unlike the pull code, the free energy
code can also handle atoms that are connected by constraints.

Potentials of mean force can also be calculated using position
restraints. With position restraints, atoms can be linked to a position
in space with a harmonic potential (see
). These positions can be
made a function of the coupling parameter :math:`\lambda`. The positions
for the A and the B states are supplied to :ref:`grompp <gmx grompp>` with
the ``-r`` and ``-rb`` options, respectively.
One could use this approach to do targeted
MD; note that we do not encourage the
use of targeted MD for proteins. A protein can be forced from one
conformation to another by using these conformations as position
restraint coordinates for state A and B. One can then slowly change
:math:`\lambda` from 0 to 1. The main drawback of this approach is that
the conformational freedom of the protein is severely limited by the
position restraints, independent of the change from state A to B. Also,
the protein is forced from state A to B in an almost straight line,
whereas the real pathway might be very different. An example of a more
fruitful application is a solid system or a liquid confined between
walls where one wants to measure the force required to change the
separation between the boundaries or walls. Because the boundaries (or
walls) already need to be fixed, the position restraints do not limit
the system in its sampling.

Removing fastest degrees of freedom
-----------------------------------

The maximum time step in MD simulations is limited by the smallest
oscillation period that can be found in the simulated system.
Bond-stretching vibrations are in their quantum-mechanical ground state
and are therefore better represented by a constraint instead of a
harmonic potential.

For the remaining degrees of freedom, the shortest oscillation period
(as measured from a simulation) is 13 fs for bond-angle vibrations
involving hydrogen atoms. Taking as a guideline that with a Verlet
(leap-frog) integration scheme a minimum of 5 numerical integration
steps should be performed per period of a harmonic oscillation in order
to integrate it with reasonable accuracy, the maximum time step will be
about 3 fs. Disregarding these very fast oscillations of period 13 fs,
the next shortest periods are around 20 fs, which will allow a maximum
time step of about 4 fs.

Removing the bond-angle degrees of freedom from hydrogen atoms can best
be done by defining them as virtual interaction
sites
instead of normal
atoms. Whereas a normal atom is connected to the molecule with bonds,
angles and dihedrals, a virtual site’s position is calculated from the
position of three nearby heavy atoms in a predefined manner (see also
sec. 
). For the hydrogens in water
and in hydroxyl, sulfhydryl, or amine groups, no degrees of freedom can
be removed, because rotational freedom should be preserved. The only
other option available to slow down these motions is to increase the
mass of the hydrogen atoms at the expense of the mass of the connected
heavy atom. This will increase the moment of inertia of the water
molecules and the hydroxyl, sulfhydryl, or amine groups, without
affecting the equilibrium properties of the system and without affecting
the dynamical properties too much. These constructions will shortly be
described in sec. 
and have previously
been described in full detail Feenstra, Hess, and Berendsen (1999).

Using both virtual sites and modified
masses
the next bottleneck is likely to be formed by the improper dihedrals
(which are used to preserve planarity or chirality of molecular groups)
and the peptide dihedrals. The peptide dihedral cannot be changed
without affecting the physical behavior of the protein. The improper
dihedrals that preserve planarity mostly deal with aromatic residues.
Bonds, angles, and dihedrals in these residues can also be replaced with
somewhat elaborate virtual site constructions.

All modifications described in this section can be performed using the
|Gromacs| topology building tool
:ref:`pdb2gmx <gmx pdb2gmx>`. Separate options
exist to increase hydrogen masses, virtualize all hydrogen atoms, or
also virtualize all aromatic residues. **Note** that when all hydrogen
atoms are virtualized, those inside the aromatic residues will be
virtualized as well, *i.e.* hydrogens in the aromatic residues are
treated differently depending on the treatment of the aromatic residues.

Parameters for the virtual site constructions for the hydrogen atoms are
inferred from the force-field parameters (*vis*. bond lengths and
angles) directly by :ref:`grompp <gmx grompp>`
while processing the topology file. The constructions for the aromatic
residues are based on the bond lengths and angles for the geometry as
described in the force fields, but these parameters are hard-coded into
:ref:`pdb2gmx <gmx pdb2gmx>` due to the complex
nature of the construction needed for a whole aromatic group.

Hydrogen bond-angle vibrations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Construction of virtual sites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fig-vsitehydro:

.. figure:: plots/dumtypes.eps

   The different types of virtual site constructions used for hydrogen
   atoms. The atoms used in the construction of the virtual site(s) are
   depicted as black circles, virtual sites as gray ones. Hydrogens are
   smaller than heavy atoms. A: fixed bond angle, note
   that here the hydrogen is not a virtual site; B: in
   the plane of three atoms, with fixed distance; C: in
   the plane of three atoms, with fixed angle and distance;
   D: construction for amine groups
   (-NH:math:`_2` or -NH:math:`_3^+`),
   see text for details.

The goal of defining hydrogen atoms as virtual sites is to remove all
high-frequency degrees of freedom from them. In some cases, not all
degrees of freedom of a hydrogen atom should be removed, *e.g.* in the
case of hydroxyl or amine groups the rotational freedom of the hydrogen
atom(s) should be preserved. Care should be taken that no unwanted
correlations are introduced by the construction of virtual sites, *e.g.*
bond-angle vibration between the constructing atoms could translate into
hydrogen bond-length vibration. Additionally, since virtual sites are by
definition massless, in order to preserve total system mass, the mass of
each hydrogen atom that is treated as virtual site should be added to
the bonded heavy atom.

Taking into account these considerations, the hydrogen atoms in a
protein naturally fall into several categories, each requiring a
different approach (see also :numref:`Fig. (%s) <fig-vsitehydro>`).

-  *hydroxyl (-OH) or sulfhydryl (-SH) hydrogen:* 
   The only internal degree of freedom in a hydroxyl group
   that can be constrained is the bending of the C-O-H
   angle. This angle is fixed by defining an additional bond of
   appropriate length, see :numref:`Fig. (%s) A<fig-vsitehydro>`.
   Doing so removes the high-frequency angle bending, but leaves the
   dihedral rotational freedom. The same goes for a sulfhydryl group.
   **Note** that in these cases the hydrogen is not treated as a virtual
   site.

-  *single amine or amide (-NH-) and aromatic hydrogens
   (-CH-):* 
   The position of these hydrogens cannot be
   constructed from a linear combination of bond vectors, because of the
   flexibility of the angle between the heavy atoms. Instead, the
   hydrogen atom is positioned at a fixed distance from the bonded heavy
   atom on a line going through the bonded heavy atom and a point on the
   line through both second bonded atoms, see
   :numref:`Fig. (%s) B<fig-vsitehydro>`.

-  *planar amine (-NH*:math:`_2`) *hydrogens:* The method
   used for the single amide hydrogen is not well suited for planar
   amine groups, because no suitable two heavy atoms can be found to
   define the direction of the hydrogen atoms. Instead, the hydrogen is
   constructed at a fixed distance from the nitrogen atom, with a fixed
   angle to the carbon atom, in the plane defined by one of the other
   heavy atoms, see :numref:`Fig. (%s) C<fig-vsitehydro>`.

-  *amine group (umbrella -NH*:math:`_2` *or
   -NH*:math:`_3^+`)* hydrogens:* Amine hydrogens with
   rotational freedom cannot be constructed as virtual sites from the
   heavy atoms they are connected to, since this would result in loss of
   the rotational freedom of the amine group. To preserve the rotational
   freedom while removing the hydrogen bond-angle degrees of freedom,
   two “dummy masses” are constructed with the same total mass, moment
   of inertia (for rotation around the C-N bond) and
   center of mass as the amine group. These dummy masses have no
   interaction with any other atom, except for the fact that they are
   connected to the carbon and to each other, resulting in a rigid
   triangle. From these three particles, the positions of the nitrogen
   and hydrogen atoms are constructed as linear combinations of the two
   carbon-mass vectors and their outer product, resulting in an amine
   group with rotational freedom intact, but without other internal
   degrees of freedom. See :numref:`Fig. (%s) D<fig-vsitehydro>`.

.. figure:: plots/dumaro.eps

   The different types of virtual site constructions used for aromatic
   residues. The atoms used in the construction of the virtual site(s)
   are depicted as black circles, virtual sites as gray ones. Hydrogens
   are smaller than heavy atoms. A: phenylalanine;
   B: tyrosine (note that the hydroxyl hydrogen is *not*
   a virtual site); C: tryptophan; D:
   histidine.

Out-of-plane vibrations in aromatic groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The planar arrangements in the side chains of the aromatic residues
lends itself perfectly to a virtual-site construction, giving a
perfectly planar group without the inherently unstable constraints that
are necessary to keep normal atoms in a plane. The basic approach is to
define three atoms or dummy masses with constraints between them to fix
the geometry and create the rest of the atoms as simple virtual sites
type (see sec. 
) from these three.
Each of the aromatic residues require a different approach:

-  *Phenylalanine:* C\ :math:`_\gamma`,
   C\ :math:`_{{\epsilon}1}`, and
   C\ :math:`_{{\epsilon}2}` are kept as normal atoms,
   but with each a mass of one third the total mass of the phenyl group.
   See :numref:`Fig. (%s) A<fig-vsitehydro>`.

-  *Tyrosine:* The ring is treated identically to the phenylalanine
   ring. Additionally, constraints are defined between
   C\ :math:`_{{\epsilon}1}`,
   C\ :math:`_{{\epsilon}2}`, and
   O\ :math:`_{\eta}`. The original improper dihedral
   angles will keep both triangles (one for the ring and one with
   O\ :math:`_{\eta}`) in a plane, but due to the larger
   moments of inertia this construction will be much more stable. The
   bond-angle in the hydroxyl group will be constrained by a constraint
   between C\ :math:`_\gamma` and
   H\ :math:`_{\eta}`. **Note** that the hydrogen is not
   treated as a virtual site. See
   :numref:`Fig. (%s) B<fig-vsitehydro>`.

-  *Tryptophan:* C\ :math:`_\beta` is kept as a normal
   atom and two dummy masses are created at the center of mass of each
   of the rings, each with a mass equal to the total mass of the
   respective ring (C\ :math:`_{{\delta}2}` and
   C\ :math:`_{{\epsilon}2}` are each counted half for
   each ring). This keeps the overall center of mass and the moment of
   inertia almost (but not quite) equal to what it was. See
   :numref:`Fig. (%s) C<fig-vsitehydro>`.

-  *Histidine:* C\ :math:`_\gamma`,
   C\ :math:`_{{\epsilon}1}` and
   N\ :math:`_{{\epsilon}2}` are kept as normal atoms,
   but with masses redistributed such that the center of mass of the
   ring is preserved. See :numref:`Fig. (%s) D<fig-vsitehydro>`.

Viscosity calculation
---------------------

The shear viscosity is a property of liquids that can be determined
easily by experiment. It is useful for parameterizing a force field
because it is a kinetic property, while most other properties which are
used for parameterization are thermodynamic. The viscosity is also an
important property, since it influences the rates of conformational
changes of molecules solvated in the liquid.

The viscosity can be calculated from an equilibrium simulation using an
Einstein relation:

.. math::

   \eta = \frac{1}{2}\frac{V}{k_B T} \lim_{t \rightarrow \infty}
   \frac{\mbox{d}}{\mbox{d} t} \left\langle 
   \left( \int_{t_0}^{{t_0}+t} P_{xz}(t') \mbox{d} t' \right)^2
   \right\rangle_{t_0}

This can be done with :ref:`gmx energy <gmx energy>`. This method converges
very slowly Hess (2002b), and as such a nanosecond simulation might not
be long enough for an accurate determination of the viscosity. The
result is very dependent on the treatment of the electrostatics. Using a
(short) cut-off results in large noise on the off-diagonal pressure
elements, which can increase the calculated viscosity by an order of
magnitude.

|Gromacs| also has a non-equilibrium method for determining the
viscosity Hess (2002b). This makes use of the fact that energy, which is
fed into system by external forces, is dissipated through viscous
friction. The generated heat is removed by coupling to a heat bath. For
a Newtonian liquid adding a small force will result in a velocity
gradient according to the following equation:

.. math:: a_x(z) + \frac{\eta}{\rho} \frac{\partial^2 v_x(z)}{\partial z^2} = 0

Here we have applied an acceleration :math:`a_x(z)` in the
:math:`x`-direction, which is a function of the :math:`z`-coordinate. In
|Gromacs| the acceleration profile is:

.. math:: a_x(z) = A \cos\left(\frac{2\pi z}{l_z}\right)

where :math:`l_z` is the height of the box. The generated velocity
profile is:

.. math:: v_x(z) = V \cos\left(\frac{2\pi z}{l_z}\right)

.. math:: V = A \frac{\rho}{\eta}\left(\frac{l_z}{2\pi}\right)^2

The viscosity can be calculated from :math:`A` and :math:`V`:

.. math:: \eta = \frac{A}{V}\rho \left(\frac{l_z}{2\pi}\right)^2
          :label: eqvisc

In the simulation :math:`V` is defined as:

.. math::

   V = \frac{\displaystyle \sum_{i=1}^N m_i v_{i,x} 2 \cos\left(\frac{2\pi z}{l_z}\right)}
            {\displaystyle \sum_{i=1}^N m_i}

The generated velocity profile is not coupled to the heat bath.
Moreover, the velocity profile is excluded from the kinetic energy. One
would like :math:`V` to be as large as possible to get good statistics.
However, the shear rate should not be so high that the system gets too
far from equilibrium. The maximum shear rate occurs where the cosine is
zero, the rate being:

.. math::

   \mbox{sh}_{\max} =  \max_z \left| \frac{\partial v_x(z)}{\partial z} \right|
   = A \frac{\rho}{\eta} \frac{l_z}{2\pi}

For a simulation with: :math:`\eta=10^{-3}`
[kgm:math:`^{-1}`\ s\ :math:`^{-1}`],
:math:`\rho=10^3`\ [kgm:math:`^{-3}`] and :math:`l_z=2\pi`\ [nm],
:math:`\mbox{sh}_{\max}=1`\ [psnm:math:`^{-1}`] :math:`A`. This shear
rate should be smaller than one over the longest correlation time in the
system. For most liquids, this will be the rotation correlation time,
which is around 10 ps. In this case, :math:`A` should be smaller than
0.1[nmps\ :math:`^{-2}`]. When the shear rate is too high, the observed
viscosity will be too low. Because :math:`V` is proportional to the
square of the box height, the optimal box is elongated in the
:math:`z`-direction. In general, a simulation length of 100 ps is enough
to obtain an accurate value for the viscosity.

The heat generated by the viscous friction is removed by coupling to a
heat bath. Because this coupling is not instantaneous the real
temperature of the liquid will be slightly lower than the observed
temperature. Berendsen derived this temperature shift Berendsen (1991),
which can be written in terms of the shear rate as:

.. math:: T_s = \frac{\eta\,\tau}{2 \rho\,C_v} \mbox{sh}_{\max}^2

where :math:`\tau` is the coupling time for the Berendsen thermostat
and :math:`C_v` is the heat capacity. Using the values of the example
above, :math:`\tau=10^{-13}` [s] and :math:`C_v=2 \cdot 10^3`\ [J
kg\ :math:`^{-1}`\ K\ :math:`^{-1}`], we get:
:math:`T_s=25`\ [Kps:math:`^{-2}`]sh\ :math:`_{\max}^2`. When we want
the shear rate to be smaller than :math:`1/10`\ [ps:math:`^{-1}`],
:math:`T_s` is smaller than 0.25[K], which is negligible.

**Note** that the system has to build up the velocity profile when
starting from an equilibrium state. This build-up time is of the order
of the correlation time of the liquid.

Two quantities are written to the energy file, along with their averages
and fluctuations: :math:`V` and :math:`1/\eta`, as obtained from
(
).

Tabulated interaction functions
-------------------------------

Cubic splines for potentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some of the inner loops of |Gromacs|, look-up tables are used for
computation of potential and forces. The tables are interpolated using a
cubic spline algorithm. There are separate tables for electrostatic,
dispersion, and repulsion interactions, but for the sake of caching
performance these have been combined into a single array. The cubic
spline interpolation for :math:`x_i \leq x < x_{i+1}` looks like this:

.. math::

   V_s(x) = A_0 + A_1 \,\epsilon + A_2 \,\epsilon^2 + A_3 \,\epsilon^3
   \label{eqn:spline}

where the table spacing :math:`h` and fraction :math:`\epsilon` are
given by:

.. math::

   \begin{aligned}
   h	&=&	x_{i+1} - x_i	\\
   \epsilon&=&	(x - x_i)/h\end{aligned}

so that :math:`0 \le \epsilon < 1`. From this, we can calculate the
derivative in order to determine the forces:

.. math::

   -V_s'(x) ~=~ 
   -\frac{{\rm d}V_s(x)}{{\rm d}\epsilon}\frac{{\rm d}\epsilon}{{\rm d}x} ~=~
   -(A_1 + 2 A_2 \,\epsilon + 3 A_3 \,\epsilon^2)/h

The four coefficients are determined from the four conditions that
:math:`V_s` and :math:`-V_s'` at both ends of each interval should match
the exact potential :math:`V` and force :math:`-V'`. This results in the
following errors for each interval:

.. math::

   \begin{aligned}
   | V_s  - V  | _{max} &=& V'''' \frac{h^4}{384} + O(h^5) \\
   | V_s' - V' | _{max} &=& V'''' \frac{h^3}{72\sqrt{3}} + O(h^4) \\
   | V_s''- V''| _{max} &=& V'''' \frac{h^2}{12}  + O(h^3)\end{aligned}

V and V’ are continuous, while V” is the first discontinuous
derivative. The number of points per nanometer is 500 and 2000 for
mixed- and double-precision versions of |Gromacs|, respectively. This
means that the errors in the potential and force will usually be smaller
than the mixed precision accuracy.

|Gromacs| stores :math:`A_0`, :math:`A_1`, :math:`A_2` and :math:`A_3`.
The force routines get a table with these four parameters and a scaling
factor :math:`s` that is equal to the number of points per nm. (**Note**
that :math:`h` is :math:`s^{-1}`). The algorithm goes a little something
like this:

#. Calculate distance vector
   (:math:`{\mbox{\boldmath ${r}$}}_{ij}`) and distance
   r\ :math:`_{ij}`

#. Multiply r\ :math:`_{ij}` by :math:`s` and truncate to an integer
   value :math:`n_0` to get a table index

#. Calculate fractional component (:math:`\epsilon` =
   :math:`s`\ r\ :math:`_{ij} - n_0`) and :math:`\epsilon^2`

#. Do the interpolation to calculate the potential :math:`V` and the
   scalar force :math:`f`

#. Calculate the vector force :math:`{\mbox{\boldmath ${F}$}}` by
   multiplying :math:`f` with
   :math:`{\mbox{\boldmath ${r}$}}_{ij}`

**Note** that table look-up is significantly *slower* than computation
of the most simple Lennard-Jones and Coulomb interaction. However, it is
much faster than the shifted Coulomb function used in conjunction with
the PPPM method. Finally, it is much easier to modify a table for the
potential (and get a graphical representation of it) than to modify the
inner loops of the MD program.

User-specified potential functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also use your own potential
functions
without editing the
|Gromacs| code. The potential function should be according to the
following equation

.. math:: V(r_{ij}) ~=~ \frac{q_i q_j}{4 \pi\epsilon_0} f(r_{ij}) + C_6 \,g(r_{ij}) + C_{12} \,h(r_{ij})

where :math:`f`, :math:`g`, and :math:`h` are user defined functions.
**Note** that if :math:`g(r)` represents a normal dispersion
interaction, :math:`g(r)` should be :math:`<` 0. C\ :math:`_6`,
C\ :math:`_{12}` and the charges are read from the topology. Also note
that combination rules are only supported for Lennard-Jones and
Buckingham, and that your tables should match the parameters in the
binary topology.

When you add the following lines in your :ref:`mdp` file:

::

    rlist           = 1.0
    coulombtype     = User
    rcoulomb        = 1.0
    vdwtype         = User
    rvdw            = 1.0

:ref:`mdrun <gmx mdrun>` will read a single non-bonded table file, or
multiple when ``energygrp-table`` is set (see below). The
name of the file(s) can be set with the :ref:`mdrun <gmx mdrun>` option
``-table``. The table file should contain seven columns of
table look-up data in the order: :math:`x`, :math:`f(x)`,
:math:`-f'(x)`, :math:`g(x)`, :math:`-g'(x)`, :math:`h(x)`,
:math:`-h'(x)`. The :math:`x` should run from 0 to :math:`r_c+1` (the
value of ``table_extension`` can be changed in the
:ref:`mdp` file). You can choose the spacing you like; for the
standard tables |Gromacs| uses a spacing of 0.002 and 0.0005 nm when you
run in mixed and double precision, respectively. In this context,
:math:`r_c` denotes the maximum of the two cut-offs ``rvdw``
and ``rcoulomb`` (see above). These variables need not be
the same (and need not be 1.0 either). Some functions used for
potentials contain a singularity at :math:`x = 0`, but since atoms are
normally not closer to each other than 0.1 nm, the function value at
:math:`x = 0` is not important. Finally, it is also possible to combine
a standard Coulomb with a modified LJ potential (or vice versa). One
then specifies *e.g.* ``coulombtype = Cut-off`` or
``coulombtype = PME``, combined with ``vdwtype
= User``. The table file must always contain the 7 columns however, and
meaningful data (i.e. not zeroes) must be entered in all columns. A
number of pre-built table files can be found in the
``GMXLIB`` directory for 6-8, 6-9, 6-10, 6-11, and 6-12
Lennard-Jones potentials combined with a normal Coulomb.

If you want to have different functional forms between different groups
of atoms, this can be set through energy groups. Different tables can be
used for non-bonded interactions between different energy groups pairs
through the :ref:`mdp` option ``energygrp-table``
(see details in the User Guide). Atoms that should interact with a
different potential should be put into different energy groups. Between
group pairs which are not listed in ``energygrp-table``, the
normal user tables will be used. This makes it easy to use a different
functional form between a few types of atoms.

Mixed Quantum-Classical simulation techniques
---------------------------------------------

In a molecular mechanics (MM) force field, the influence of electrons is
expressed by empirical parameters that are assigned on the basis of
experimental data, or on the basis of results from high-level quantum
chemistry calculations. These are valid for the ground state of a given
covalent structure, and the MM approximation is usually sufficiently
accurate for ground-state processes in which the overall connectivity
between the atoms in the system remains unchanged. However, for
processes in which the connectivity does change, such as chemical
reactions, or processes that involve multiple electronic states, such as
photochemical conversions, electrons can no longer be ignored, and a
quantum mechanical description is required for at least those parts of
the system in which the reaction takes place.

One approach to the simulation of chemical reactions in solution, or in
enzymes, is to use a combination of quantum mechanics (QM) and molecular
mechanics (MM). The reacting parts of the system are treated quantum
mechanically, with the remainder being modeled using the force field.
The current version of |Gromacs| provides interfaces to several popular
Quantum Chemistry packages (MOPAC Dewar (1983), GAMESS-UK Guest et al.
(1987), Gaussian Frisch et al. (n.d.) and CPMD Car and Parrinello
(1985)).

|Gromacs| interactions between the two subsystems are either handled as
described by Field et al. Field, Bash, and Karplus (1990) or within
the ONIOM approach by Morokuma and coworkers Maseras and Morokuma (1995;
Svensson et al. 1996).

Overview
^^^^^^^^

Two approaches for describing the interactions between the QM and MM
subsystems are supported in this version:

#. **Electronic Embedding** The electrostatic interactions between the
   electrons of the QM region and the MM atoms and between the QM nuclei
   and the MM atoms are included in the Hamiltonian for the QM
   subsystem:

   .. math::

      H^{QM/MM} =
      H^{QM}_e-\sum_i^n\sum_J^M\frac{e^2Q_J}{4\pi\epsilon_0r_{iJ}}+\sum_A^N\sum_J^M\frac{e^2Z_AQ_J}{e\pi\epsilon_0R_{AJ}},

#  where :math:`n` and :math:`N` are the number of electrons and nuclei
   in the QM region, respectively, and :math:`M` is the number of
   charged MM atoms. The first term on the right hand side is the
   original electronic Hamiltonian of an isolated QM system. The first
   of the double sums is the total electrostatic interaction between the
   QM electrons and the MM atoms. The total electrostatic interaction of
   the QM nuclei with the MM atoms is given by the second double sum.
   Bonded interactions between QM and MM atoms are described at the MM
   level by the appropriate force-field terms. Chemical bonds that
   connect the two subsystems are capped by a hydrogen atom to complete
   the valence of the QM region. The force on this atom, which is
   present in the QM region only, is distributed over the two atoms of
   the bond. The cap atom is usually referred to as a link atom.

#. **ONIOM** In the ONIOM approach, the energy and gradients are first
   evaluated for the isolated QM subsystem at the desired level of *ab
   initio* theory. Subsequently, the energy and gradients of the total
   system, including the QM region, are computed using the molecular
   mechanics force field and added to the energy and gradients
   calculated for the isolated QM subsystem. Finally, in order to
   correct for counting the interactions inside the QM region twice, a
   molecular mechanics calculation is performed on the isolated QM
   subsystem and the energy and gradients are subtracted. This leads to
   the following expression for the total QM/MM energy (and gradients
   likewise):

   .. math::

      E_{tot} = E_{I}^{QM}
      +E_{I+II}^{MM}-E_{I}^{MM},

#  where the subscripts I and II refer to the QM and MM subsystems,
   respectively. The superscripts indicate at what level of theory the
   energies are computed. The ONIOM scheme has the advantage that it is
   not restricted to a two-layer QM/MM description, but can easily
   handle more than two layers, with each layer described at a different
   level of theory.

Usage
^^^^^

To make use of the QM/MM functionality in |Gromacs|, one needs to:

#. introduce link atoms at the QM/MM boundary, if needed;

#. specify which atoms are to be treated at a QM level;

#. specify the QM level, basis set, type of QM/MM interface and so on.

Adding link atoms
^^^^^^^^^^^^^^^^^

At the bond that connects the QM and MM subsystems, a link atoms is
introduced. In |Gromacs| the link atom has special atomtype, called LA.
This atomtype is treated as a hydrogen atom in the QM calculation, and
as a virtual site in the force-field calculation. The link atoms, if
any, are part of the system, but have no interaction with any other
atom, except that the QM force working on it is distributed over the two
atoms of the bond. In the topology, the link atom (LA), therefore, is
defined as a virtual site atom:

::

    [ virtual_sites2 ]
    LA QMatom MMatom 1 0.65

See sec. 
for more details on how virtual
sites are treated. The link atom is replaced at every step of the
simulation.

In addition, the bond itself is replaced by a constraint:

::

    [ constraints ]
    QMatom MMatom 2 0.153

**Note** that, because in our system the QM/MM bond is a carbon-carbon
bond (0.153 nm), we use a constraint length of 0.153 nm, and dummy
position of 0.65. The latter is the ratio between the ideal C-H bond
length and the ideal C-C bond length. With this ratio, the link atom is
always 0.1 nm away from the ``QMatom``, consistent with the
carbon-hydrogen bond length. If the QM and MM subsystems are connected
by a different kind of bond, a different constraint and a different
dummy position, appropriate for that bond type, are required.

Specifying the QM atoms
^^^^^^^^^^^^^^^^^^^^^^^

Atoms that should be treated at a QM level of theory, including the link
atoms, are added to the index file. In addition, the chemical bonds
between the atoms in the QM region are to be defined as connect bonds
(bond type 5) in the topology file:

::

    [ bonds ]
    QMatom1 QMatom2 5
    QMatom2 QMatom3 5

Specifying the QM/MM simulation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the :ref:`mdp` file, the following parameters control a
QM/MM simulation.

``QMMM = no``
    | If this is set to ``yes``, a QM/MM simulation is
      requested. Several groups of atoms can be described at different
      QM levels separately. These are specified in the QMMM-grps field
      separated by spaces. The level of *ab initio* theory at which the
      groups are described is specified by ``QMmethod`` and
      ``QMbasis`` Fields. Describing the groups at different
      levels of theory is only possible with the ONIOM QM/MM scheme,
      specified by ``QMMMscheme``.

``QMMM-grps =``
    | groups to be described at the QM level

``QMMMscheme = normal``
    | Options are ``normal`` and ``ONIOM``. This
      selects the QM/MM interface. ``normal`` implies that
      the QM subsystem is electronically embedded in the MM subsystem.
      There can only be one ``QMMM-grps`` that is modeled at
      the ``QMmethod`` and ``QMbasis`` level of
      * ab initio* theory. The rest of the system is described at the MM
      level. The QM and MM subsystems interact as follows: MM point
      charges are included in the QM one-electron Hamiltonian and all
      Lennard-Jones interactions are described at the MM level. If
      ``ONIOM`` is selected, the interaction between the
      subsystem is described using the ONIOM method by Morokuma and
      co-workers. There can be more than one QMMM-grps each modeled at a
      different level of QM theory (QMmethod and QMbasis).

``QMmethod =``
    | Method used to compute the energy and gradients on the QM atoms.
      Available methods are AM1, PM3, RHF, UHF, DFT, B3LYP, MP2, CASSCF,
      MMVB and CPMD. For CASSCF, the number of electrons and orbitals
      included in the active space is specified by
      ``CASelectrons`` and ``CASorbitals``. For
      CPMD, the plane-wave cut-off is specified by the
      ``planewavecutoff`` keyword.

``QMbasis =``
    | Gaussian basis set used to expand the electronic wave-function.
      Only Gaussian basis sets are currently available, i.e. STO-3G,
      3-21G, 3-21G\*, 3-21+G\*, 6-21G, 6-31G, 6-31G\*, 6-31+G\*, and
      6-311G. For CPMD, which uses plane wave expansion rather than
      atom-centered basis functions, the ``planewavecutoff``
      keyword controls the plane wave expansion.

``QMcharge =``
    | The total charge in *e* of the ``QMMM-grps``. In case
      there are more than one ``QMMM-grps``, the total
      charge of each ONIOM layer needs to be specified separately.

``QMmult =``
    | The multiplicity of the ``QMMM-grps``. In case there
      are more than one ``QMMM-grps``, the multiplicity of
      each ONIOM layer needs to be specified separately.

``CASorbitals =``
    | The number of orbitals to be included in the active space when
      doing a CASSCF computation.

``CASelectrons =``
    | The number of electrons to be included in the active space when
      doing a CASSCF computation.

``SH = no``
    | If this is set to yes, a QM/MM MD simulation on the excited
      state-potential energy surface and enforce a diabatic hop to the
      ground-state when the system hits the conical intersection
      hyperline in the course the simulation. This option only works in
      combination with the CASSCF method.

Output
^^^^^^

The energies and gradients computed in the QM calculation are added to
those computed by |Gromacs|. In the :ref:`edr` file there is a
section for the total QM energy.

Future developments
^^^^^^^^^^^^^^^^^^^

Several features are currently under development to increase the
accuracy of the QM/MM interface. One useful feature is the use of
delocalized MM charges in the QM computations. The most important
benefit of using such smeared-out charges is that the Coulombic
potential has a finite value at interatomic distances. In the point
charge representation, the partially-charged MM atoms close to the QM
region tend to “over-polarize” the QM system, which leads to artifacts
in the calculation.

What is needed as well is a transition state optimizer.

Using VMD plug-ins for trajectory file I/O
------------------------------------------

|Gromacs|
tools are able to use the plug-ins found in an existing installation of
`VMD <http://www.ks.uiuc.edu/Research/vmd>`__ in order to read and write
trajectory files in formats that are not native to |Gromacs|. You will be
able to supply an AMBER DCD-format trajectory filename directly to
|Gromacs| tools, for example.

This requires a VMD installation not older than version 1.8, that your
system provides the dlopen function so that programs can determine at
run time what plug-ins exist, and that you build shared libraries when
building |Gromacs|. CMake will find the vmd executable in your path, and
from it, or the environment variable ``VMDDIR`` at
configuration or run time, locate the plug-ins. Alternatively, the
``VMD_PLUGIN_PATH`` can be used at run time to specify a
path where these plug-ins can be found. Note that these plug-ins are in
a binary format, and that format must match the architecture of the
machine attempting to use them.

Interactive Molecular Dynamics
------------------------------

|Gromacs| supports the interactive molecular dynamics (IMD) protocol as
implemented by `VMD <http://www.ks.uiuc.edu/Research/vmd>`__ to control
a running simulation in NAMD. IMD allows to monitor a running |Gromacs|
simulation from a VMD client. In addition, the user can interact with
the simulation by pulling on atoms, residues or fragments with a mouse
or a force-feedback device. Additional information about the |Gromacs|
implementation and an exemplary |Gromacs| IMD system can be found `on this
homepage <http://www.mpibpc.mpg.de/grubmueller/interactivemd>`__.

Simulation input preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The |Gromacs| implementation allows transmission and interaction with a
part of the running simulation only, e.g. in cases where no water
molecules should be transmitted or pulled. The group is specified via
the :ref:`mdp` option ``IMD-group``. When
``IMD-group`` is empty, the IMD protocol is disabled and
cannot be enabled via the switches in :ref:`mdrun <gmx mdrun>`. To interact
with the entire system, ``IMD-group`` can be set to
``System``. When using :ref:`grompp <gmx grompp>`, a
:ref:`gro` file to be used as VMD input is written out
(``-imd`` switch of :ref:`grompp <gmx grompp>`).

Starting the simulation
^^^^^^^^^^^^^^^^^^^^^^^

Communication between VMD and |Gromacs| is achieved via TCP sockets and
thus enables controlling an :ref:`mdrun <gmx mdrun>` running locally or on
a remote cluster. The port for the connection can be specified with the
``-imdport`` switch of :ref:`mdrun <gmx mdrun>`, 8888 is the
default. If a port number of 0 or smaller is provided, |Gromacs|
automatically assigns a free port to use with IMD.

Every :math:`N` steps, the :ref:`mdrun <gmx mdrun>` client receives the
applied forces from VMD and sends the new positions to the client. VMD
permits increasing or decreasing the communication frequency
interactively. By default, the simulation starts and runs even if no IMD
client is connected. This behavior is changed by the
``-imdwait`` switch of :ref:`mdrun <gmx mdrun>`. After startup
and whenever the client has disconnected, the integration stops until
reconnection of the client. When the ``-imdterm`` switch is
used, the simulation can be terminated by pressing the stop button in
VMD. This is disabled by default. Finally, to allow interacting with the
simulation (i.e. pulling from VMD) the ``-imdpull`` switch
has to be used. Therefore, a simulation can only be monitored but not
influenced from the VMD client when none of ``-imdwait``,
``-imdterm`` or ``-imdpull`` are set. However,
since the IMD protocol requires no authentication, it is not advisable
to run simulations on a host directly reachable from an insecure
environment. Secure shell forwarding of TCP can be used to connect to
running simulations not directly reachable from the interacting host.
Note that the IMD command line switches of :ref:`mdrun <gmx mdrun>` are
hidden by default and show up in the help text only with
:ref:`gmx mdrun` ``-h -hidden``.

Connecting from VMD
^^^^^^^^^^^^^^^^^^^

In VMD, first the structure corresponding to the IMD group has to be
loaded (*File* :math:`\rightarrow` *New Molecule*). Then the IMD
connection window has to be used (*Extensions* :math:`\rightarrow`
*Simulation* :math:`\rightarrow` *IMD Connect (NAMD)*). In the IMD
connection window, hostname and port have to be specified and followed
by pressing *Connect*. *Detach Sim* allows disconnecting without
terminating the simulation, while *Stop Sim* ends the simulation on the
next neighbor searching step (if allowed by ``-imdterm``).

The timestep transfer rate allows adjusting the communication frequency
between simulation and IMD client. Setting the keep rate loads every
:math:`N^\mathrm{th}` frame into VMD instead of discarding them when a
new one is received. The displayed energies are in SI units in contrast
to energies displayed from NAMD simulations.

Embedding proteins into the membranes
-------------------------------------

|Gromacs| is capable of inserting the protein into pre-equilibrated lipid
bilayers with minimal perturbation of the lipids using the method, which
was initially described as a ProtSqueeze technique,Yesylevskyy (2007)
and later implemented as g_membed tool. Wolf et al. (2010) Currently the
functionality of g_membed is available in mdrun as described in the
user guide.

This method works by first artificially shrinking the protein in the
:math:`xy`-plane, then it removes lipids that overlap with that much
smaller core. Then the protein atoms are gradually resized back to their
initial configuration, using normal dynamics for the rest of the system,
so the lipids adapt to the protein. Further lipids are removed as
required.
