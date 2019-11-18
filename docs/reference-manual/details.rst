Some implementation details
===========================

In this chapter we will present some implementation details. This is far
from complete, but we deemed it necessary to clarify some things that
would otherwise be hard to understand.

Single Sum Virial in |Gromacs|
------------------------------

The virial :math:`\Xi` can be written in full tensor form as:

.. math:: \Xi~=~-\frac{1}{2}~\sum_{i < j}^N~\mathbf{r}_{ij}\otimes\mathbf{F}_{ij}
          :label: eqnvirialfulltensorform

where :math:`\otimes` denotes the *direct product* of two vectors. [1]_
When this is computed in the inner loop of an MD program 9
multiplications and 9 additions are needed. [2]_

Here it is shown how it is possible to extract the virial calculation
from the inner loop \ :ref:`177 <refBekker93b>`.

Virial
~~~~~~

In a system with periodic boundary conditions, the periodicity must be
taken into account for the virial:

.. math:: \Xi~=~-\frac{1}{2}~\sum_{i < j}^{N}~\mathbf{r}_{ij}^n\otimes\mathbf{F}_{ij}
          :label: eqnvirialperiodicity

where :math:`\mathbf{r}_{ij}^n` denotes the distance
vector of the *nearest image* of atom :math:`i` from atom :math:`j`. In
this definition we add a *shift vector* :math:`\delta_i` to the position
vector :math:`\mathbf{r}_i` of atom :math:`i`. The
difference vector :math:`\mathbf{r}_{ij}^n` is thus equal
to:

.. math:: \mathbf{r}_{ij}^n~=~\mathbf{r}_i+\delta_i-\mathbf{r}_j
          :label: eqnvirialdiffvector

or in shorthand:

.. math:: \mathbf{r}_{ij}^n~=~\mathbf{r}_i^n-\mathbf{r}_j
          :label: eqnvirialdiffvecshort

In a triclinic system, there are 27 possible images of :math:`i`; when
a truncated octahedron is used, there are 15 possible images.

Virial from non-bonded forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here the derivation for the single sum virial in the *non-bonded force*
routine is given. There are a couple of considerations that are special
to |Gromacs| that we take into account:

-  When calculating short-range interactions, we apply the *minimum
   image convention* and only consider the closest image of each
   neighbor - and in particular we never allow interactions between a
   particle and any of its periodic images. For all the equations below,
   this means :math:`i \neq j`.

-  In general, either the :math:`i` or :math:`j` particle might be
   shifted to a neighbor cell to get the closest interaction (shift
   :math:`\delta_{ij}`). However, with minimum image convention there
   can be at most 27 different shifts for particles in the central cell,
   and for typical (very short-ranged) biomolecular interactions there
   are typically only a few different shifts involved for each particle,
   not to mention that each interaction can only be present for one
   shift.

-  For the |Gromacs| nonbonded interactions we use this to split the
   neighborlist of each :math:`i` particle into multiple separate lists,
   where each list has a constant shift :math:`\delta_i` for the
   :math:`i` partlcle. We can represent this as a sum over shifts (for
   which we use index :math:`s`), with the constraint that each particle
   interaction can only contribute to one of the terms in this sum, and
   the shift is no longer dependent on the :math:`j` particles. For any
   sum that does not contain complex dependence on :math:`s`, this means
   the sum trivially reduces to just the sum over :math:`i` and/or
   :math:`j`.

-  To simplify some of the sums, we replace sums over :math:`j<i` with
   double sums over all particles (remember, :math:`i \neq j`) and
   divide by 2.

Starting from the above definition of the virial, we then get

.. math:: \begin{aligned}
          \Xi
          &~=~&-{\frac{1}{2}}~\sum_{i < j}^{N}~{\mathbf r}^n_{ij} \otimes {\mathbf F}_{ij} \nonumber \\
          &~=~&-{\frac{1}{2}}~\sum_{i < j}^{N}~\left( {\mathbf r}_i + \delta_{ij} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij} \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{j=1}^{N}~\left( {\mathbf r}_i + \delta_{ij} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij} \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij,s} \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{i=}^{N}~\sum_{s}~\sum_{j=1}^{N}~\left( \left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} -{\mathbf r}_j \otimes {\mathbf F}_{ij,s} \right) \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^N ~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^{N} {\mathbf r}_j \otimes {\mathbf F}_{ij,s} \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^N ~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{i=1}^{N}~\sum_{j=1}^{N} {\mathbf r}_j \otimes {\mathbf F}_{ij} \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes ~\sum_{j=1}^N {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{j=1}^N {\mathbf r}_j \otimes \sum_{i=1}^{N} {\mathbf F}_{ij} \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes ~\sum_{j=1}^N {\mathbf F}_{ij,s} - {\frac{1}{4}}\sum_{j=1}^N {\mathbf r}_j \otimes \sum_{i=1}^{N} {\mathbf F}_{ji} \nonumber \\
          &~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{i,s} - {\frac{1}{4}}\sum_{j=1}^N~{\mathbf r}_j \otimes {\mathbf F}_{j}  \nonumber \\
          &~=~&-{\frac{1}{4}}~\left(\sum_{i=1}^{N}~{\mathbf r}_i  \otimes {\mathbf F}_{i} + \sum_{j=1}^N~{\mathbf r}_j \otimes {\mathbf F}_{j} \right) - {\frac{1}{4}}\sum_{s}~\sum_{i=1}^{N} \delta_{i,s} \otimes {\mathbf F}_{i,s}  \nonumber \\
          &~=~&-{\frac{1}{2}}\sum_{i=1}^{N}~{\mathbf r}_i \otimes {\mathbf F}_{i} -{\frac{1}{4}}\sum_{s}~\sum_{i=1}^{N}~\delta_{i,s} \otimes {\mathbf F}_{i,s} \nonumber \\
          &~=~&-{\frac{1}{2}}\sum_{i=1}^{N}~{\mathbf r}_i \otimes {\mathbf F}_{i} -{\frac{1}{4}}\sum_{s}~\delta_{s} \otimes {\mathbf F}_{s} \nonumber \\
          &~=~&\Xi_0 + \Xi_1\end{aligned}
          :label: eqnviriallong

In the second-last stage, we have used the property that each shift
vector itself does not depend on the coordinates of particle :math:`i`,
so it is possible to sum up all forces corresponding to each shift
vector (in the nonbonded kernels), and then just use a sum over the
different shift vectors outside the kernels. We have also used

.. math:: \begin{aligned}
          \mathbf{F}_i&~=~&\sum_{j=1}^N~\mathbf{F}_{ij}					\\
          \mathbf{F}_j&~=~&\sum_{i=1}^N~\mathbf{F}_{ji}\end{aligned}
          :label: eqnvirialtotalforce

which is the total force on :math:`i` with respect to :math:`j`.
Because we use Newton’s Third Law:

.. math:: \mathbf{F}_{ij}~=~-\mathbf{F}_{ji}
          :label: eqnnewtonsthird

we must, in the implementation, double the term containing the shift
:math:`\delta_i`. Similarly, in a few places we have summed the
shift-dependent force over all shifts to come up with the total force
per interaction or particle.

This separates the total virial :math:`\Xi` into a component
:math:`\Xi_0` that is a single sum over particles, and a second
component :math:`\Xi_1` that describes the influence of the particle
shifts, and that is only a sum over the different shift vectors.

The intra-molecular shift (mol-shift)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the bonded forces and SHAKE it is possible to make a *mol-shift*
list, in which the periodicity is stored. We simple have an array mshift
in which for each atom an index in the shiftvec array is stored.

The algorithm to generate such a list can be derived from graph theory,
considering each particle in a molecule as a bead in a graph, the bonds
as edges.

#. Represent the bonds and atoms as bidirectional graph

#. Make all atoms white

#. Make one of the white atoms black (atom :math:`i`) and put it in the
   central box

#. Make all of the neighbors of :math:`i` that are currently white, gray

#. Pick one of the gray atoms (atom :math:`j`), give it the correct
   periodicity with respect to any of its black neighbors and make it
   black

#. Make all of the neighbors of :math:`j` that are currently white, gray

#. If any gray atom remains, go to [5]

#. If any white atom remains, go to [3]

Using this algorithm we can

-  optimize the bonded force calculation as well as SHAKE

-  calculate the virial from the bonded forces in the single sum method
   again

Find a representation of the bonds as a bidirectional graph.

Virial from Covalent Bonds
~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the covalent bond force gives a contribution to the virial, we
have:

.. math:: \begin{aligned}
          b	&~=~&	\|\mathbf{r}_{ij}^n\|					\\
          V_b	&~=~&	\frac{1}{2} k_b(b-b_0)^2				\\
          \mathbf{F}_i	&~=~&	-\nabla V_b					\\
          	&~=~&	k_b(b-b_0)\frac{\mathbf{r}_{ij}^n}{b}			\\
          \mathbf{F}_j	&~=~&	-\mathbf{F}_i\end{aligned}
          :label: eqncovbondvirial

The virial contribution from the bonds then is:

.. math:: \begin{aligned}
          \Xi_b	&~=~&	-\frac{1}{2}(\mathbf{r}_i^n\otimes\mathbf{F}_i~+~\mathbf{r}_j\otimes\mathbf{F}_j)	\\
          &~=~&	-\frac{1}{2}\mathbf{r}_{ij}^n\otimes\mathbf{F}_i\end{aligned}
          :label: eqncovbondvirialcontri

Virial from SHAKE
~~~~~~~~~~~~~~~~~

An important contribution to the virial comes from shake. Satisfying the
constraints a force **G** that is exerted on the particles “shaken.” If
this force does not come out of the algorithm (as in standard SHAKE) it
can be calculated afterward (when using *leap-frog*) by:

.. math:: \begin{aligned}
          \Delta\mathbf{r}_i&~=~&{\mathbf{r}_i}(t+{\Delta t})-
          [\mathbf{r}_i(t)+{\bf v}_i(t-\frac{{\Delta t}}{2}){\Delta t}+\frac{\mathbf{F}_i}{m_i}{\Delta t}^2]	\\
          {\bf G}_i&~=~&\frac{m_i{\Delta}{\mathbf{r}_i}}{{\Delta t}^2i}\end{aligned}
          :label: eqnshakevirial

This does not help us in the general case. Only when no periodicity is
needed (like in rigid water) this can be used, otherwise we must add the
virial calculation in the inner loop of SHAKE.

When it *is* applicable the virial can be calculated in the single sum
way:

.. math:: \Xi~=~-\frac{1}{2}\sum_i^{N_c}~\mathbf{r}_i\otimes\mathbf{F}_i
          :label: eqnshakevirialsinglesum

where :math:`N_c` is the number of constrained atoms.

Optimizations
-------------

Here we describe some of the algorithmic optimizations used in |Gromacs|,
apart from parallelism.

.. _waterloops:

Inner Loops for Water
~~~~~~~~~~~~~~~~~~~~~

|Gromacs| uses special inner loops to calculate non-bonded interactions
for water molecules with other atoms, and yet another set of loops for
interactions between pairs of water molecules. There highly optimized
loops for two types of water models. For three site models similar to
SPC \ :ref:`80 <refBerendsen81>`, *i.e.*:

#. There are three atoms in the molecule.

#. The whole molecule is a single charge group.

#. The first atom has Lennard-Jones (sec. :ref:`lj`) and Coulomb
   (sec. :ref:`coul`) interactions.

#. Atoms two and three have only Coulomb interactions, and equal
   charges.

These loops also works for the SPC/E \ :ref:`178 <refBerendsen87>` and
TIP3P \ :ref:`128 <refJorgensen83>` water models. And for four site water
models similar to TIP4P \ :ref:`128 <refJorgensen83>`:

#. There are four atoms in the molecule.

#. The whole molecule is a single charge group.

#. The first atom has only Lennard-Jones (sec. :ref:`lj`) interactions.

#. Atoms two and three have only Coulomb (sec. :ref:`coul`) interactions,
   and equal charges.

#. Atom four has only Coulomb interactions.

The benefit of these implementations is that there are more
floating-point operations in a single loop, which implies that some
compilers can schedule the code better. However, it turns out that even
some of the most advanced compilers have problems with scheduling,
implying that manual tweaking is necessary to get optimum performance.
This may include common-sub-expression elimination, or moving code
around.

.. raw:: latex

    \clearpage


.. [1]
   Note that some derivations, an alternative notation
   :math:`\xi_{\mathrm{alt}} = v_{\xi} = p_{\xi}/Q` is used.

.. [2]
   The calculation of Lennard-Jones and Coulomb forces is about 50
   floating point operations.
