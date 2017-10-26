Some implementation details
===========================

In this chapter we will present some implementation details. This is far
from complete, but we deemed it necessary to clarify some things that
would otherwise be hard to understand.

Single Sum Virial in GROMACS
----------------------------

The virial:raw-latex:`\index{virial}` :math:`\Xi` can be written in full
tensor form as:

.. math:: \Xi~=~-\half~\sum_{i < j}^N~\ve{r}_ij\otimes\ve{F}_ij

 where :math:`\otimes` denotes the *direct product* of two
vectors.:raw-latex:`\mbox{\boldmath ${f}$}`otnote
:math:`({\bf u}\otimes{\bf v})^{{\alpha\beta}}~=~{\bf u}_{{\alpha}}{\bf v}_{{\beta}}`
When this is computed in the inner loop of an MD program 9
multiplications and 9 additions are needed. [3]_

Here it is shown how it is possible to extract the virial calculation
from the inner loop Bekker, Berendsen, Dijkstra, Achterop, Drunen, et
al. (1993a).

Virial
~~~~~~

In a system with periodic boundary
conditions:raw-latex:`\index{periodic boundary 
conditions}`, the periodicity must be taken into account for the virial:

.. math:: \Xi~=~-\half~\sum_{i < j}^{N}~\ve{r}_{ij}^n\otimes\ve{F}_ij

 where :math:`{\mbox{\boldmath ${r}$}}_{ij}^n` denotes the distance
vector of the *nearest image* of atom :math:`i` from atom :math:`j`. In
this definition we add a *shift vector* :math:`\delta_i` to the position
vector :math:`{\mbox{\boldmath ${r}$}}_i` of atom :math:`i`. The
difference vector :math:`{\mbox{\boldmath ${r}$}}_{ij}^n` is thus equal
to:

.. math:: {\mbox{\boldmath ${r}$}}_{ij}^n~=~{\mbox{\boldmath ${r}$}}_i+\delta_i-{\mbox{\boldmath ${r}$}}_j

 or in shorthand:

.. math:: {\mbox{\boldmath ${r}$}}_{ij}^n~=~{\mbox{\boldmath ${r}$}}_i^n-{\mbox{\boldmath ${r}$}}_j

 In a triclinic system, there are 27 possible images of :math:`i`; when
a truncated octahedron is used, there are 15 possible images.

Virial from non-bonded forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here the derivation for the single sum virial in the *non-bonded force*
routine is given. :math:`i \neq j` in all formulae below.

.. math::

   \begin{aligned}
   \Xi	
   &~=~&-\half~\sum_{i < j}^{N}~\ve{r}_{ij}^n\otimes\ve{F}_ij				\\
   &~=~&-\qrt\sum_{i=1}^N~\sum_{j=1}^N ~(\ve{r}_i+\di-\ve{r}_j)\otimes\ve{F}_ij	\\
   &~=~&-\qrt\sum_{i=1}^N~\sum_{j=1}^N ~(\ve{r}_i+\di)\otimes\ve{F}_ij-\ve{r}_j\otimes\ve{F}_{ij}	\\
   &~=~&-\qrt\left(\sum_{i=1}^N~\sum_{j=1}^N ~(\ve{r}_i+\di)\otimes\ve{F}_ij~-~\sum_{i=1}^N~\sum_{j=1}^N ~\ve{r}_j\otimes\ve{F}_{ij}\right)	\\
   &~=~&-\qrt\left(\sum_{i=1}^N~(\ve{r}_i+\di)\otimes\sum_{j=1}^N~\ve{F}_ij~-~\sum_{j=1}^N ~\ve{r}_j\otimes\sum_{i=1}^N~\ve{F}_{ij}\right)	\\
   &~=~&-\qrt\left(\sum_{i=1}^N~(\ve{r}_i+\di)\otimes\ve{F}_i~+~\sum_{j=1}^N ~\ve{r}_j\otimes\ve{F}_j\right)	\\
   &~=~&-\qrt\left(2~\sum_{i=1}^N~\ve{r}_i\otimes\ve{F}_i+\sum_{i=1}^N~\di\otimes\Fvi\right)\end{aligned}

In these formulae we introduced:

.. math::

   \begin{aligned}
   \ve{F}_i&~=~&\sum_{j=1}^N~\ve{F}_{ij}					\\
   \ve{F}_j&~=~&\sum_{i=1}^N~\ve{F}_{ji}\end{aligned}

 which is the total force on :math:`i` with respect to :math:`j`.
Because we use Newton’s Third Law:

.. math:: \ve{F}_ij~=~-\ve{F}_ji

 we must, in the implementation, double the term containing the shift
:math:`\delta_i`.

The intra-molecular shift (mol-shift)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the bonded forces and SHAKE it is possible to make a *mol-shift*
list, in which the periodicity is stored. We simple have an array
:raw-latex:`\tt`mshift in which for each atom an index in the
:raw-latex:`\tt`shiftvec array is stored.

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

.. math::

   \begin{aligned}
   b	&~=~&	\|\ve{r}_{ij}^n\|					\\
   V_b	&~=~&	\half k_b(b-b_0)^2				\\
   \ve{F}_i	&~=~&	-\nabla V_b					\\
   	&~=~&	k_b(b-b_0)\frac{\ve{r}_{ij}^n}{b}			\\
   \ve{F}_j	&~=~&	-\ve{F}_i\end{aligned}

 The virial contribution from the bonds then is:

.. math::

   \begin{aligned}
   \Xi_b	&~=~&	-\half(\ve{r}_i^n\otimes\ve{F}_i~+~\ve{r}_j\otimes\ve{F}_j)	\\
   	&~=~&	-\half\ve{r}_{ij}^n\otimes\ve{F}_i\end{aligned}

Virial from SHAKE
~~~~~~~~~~~~~~~~~

An important contribution to the virial comes from shake. Satisfying the
constraints a force **G** that is exerted on the particles “shaken.” If
this force does not come out of the algorithm (as in standard SHAKE) it
can be calculated afterward (when using *leap-frog*) by:

.. math::

   \begin{aligned}
   \Delta\ve{r}_i&~=~&\rvi(t+\Dt)-
   [\ve{r}_i(t)+{\bf v}_i(t-\frac{\Dt}{2})\Dt+\frac{\ve{F}_i}{m_i}\Dt^2]	\\
   {\bf G}_i&~=~&\frac{m_i\Delta\ve{r}_i}{\Dt^2}\end{aligned}

 This does not help us in the general case. Only when no periodicity is
needed (like in rigid water) this can be used, otherwise we must add the
virial calculation in the inner loop of SHAKE.

When it *is* applicable the virial can be calculated in the single sum
way:

.. math:: \Xi~=~-\half\sum_i^{N_c}~\ve{r}_i\otimes\ve{F}_i

 where :math:`N_c` is the number of constrained atoms.

Optimizations
-------------

Here we describe some of the algorithmic optimizations used in GROMACS,
apart from parallelism.

Inner Loops for Water
~~~~~~~~~~~~~~~~~~~~~

GROMACS uses special inner loops to calculate non-bonded interactions
for water molecules with other atoms, and yet another set of loops for
interactions between pairs of water molecules. There highly optimized
loops for two types of water models. For three site models similar to
SPC Berendsen et al. (1981), *i.e.*:

#. There are three atoms in the molecule.

#. The whole molecule is a single charge group.

#. The first atom has Lennard-Jones (sec. :raw-latex:`\ref{sec:lj}`) and
   Coulomb (sec. :raw-latex:`\ref{sec:coul}`) interactions.

#. Atoms two and three have only Coulomb interactions, and equal
   charges.

These loops also works for the SPC/E Berendsen, Grigera, and Straatsma
(1987) and TIP3P Jorgensen et al. (1983) water models. And for four site
water models similar to TIP4P Jorgensen et al. (1983):

#. There are four atoms in the molecule.

#. The whole molecule is a single charge group.

#. The first atom has only Lennard-Jones
   (sec. :raw-latex:`\ref{sec:lj}`) interactions.

#. Atoms two and three have only Coulomb
   (sec. :raw-latex:`\ref{sec:coul}`) interactions, and equal charges.

#. Atom four has only Coulomb interactions.

The benefit of these implementations is that there are more
floating-point operations in a single loop, which implies that some
compilers can schedule the code better. However, it turns out that even
some of the most advanced compilers have problems with scheduling,
implying that manual tweaking is necessary to get optimum
performance:raw-latex:`\index{performance}`. This may include
common-sub-expression elimination, or moving code around.
