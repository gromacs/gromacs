Design goals and motivation for the data format of bonded forces in NB-LIB
--------------------------------------------------------------------------


The current format for listed forces in |Gromacs| looks like this:

.. code:: cpp

   struct InteractionDefinitions
   {
       std::vector<t_iparams> iparams;
       std::array<std::vector<int>, F_NRE> il;
   };

The format covers all interaction types, i.e. \ ``t_iparams`` is a union
type which can hold the parameters of any type.
The other member called ``il`` contains the
indices for each interaction type, where ``F_NRE`` is the number of
interaction types that |Gromacs| supports. More precisely, each
member of ``il``, a ``std::vector<int>``, is a flattened list of all
interactions for a given interaction type. The vector contains ``N+1`` integer indices
for each interaction, where ``N`` is the number of particles that are
involved in the interaction. An additional index is needed to retrieve
the correct parameters in ``iparams``, hence the total number of indices sums up
to ``N+1`` per interaction.

The big advantage of storing all types in a union data type is (was),
that it allows looping over all types with a simple for-loop.
In pre C++11 and perhaps even pre C++14 times, looping over different
types was a big hassle and the union data type approach likely was the
only practicable solution. One downside of this approach, however, is
that with just a single (union) type, one can't leverage the compiler's
type system, most importantly static branching, for example with overload resolution.
As a consequence, only dynamic branching with ``if`` statements remains.

Consider, for instance, the implementation of the top-level
``calc_listed(const InteractionDefinitions& idef, ...)`` in |Gromacs|, which in its essence,
looks like this:

.. code:: cpp

   void calc_listed(const InteractionDefinitions& idef, ...)
   {
       // manage timing and multi-threading 

       for (int ftype = 0; ftype < F_NRE; ++type)
       {
           // branch out and descend stack for 2 intermediate functions based on
           // the type of interaction that ftype corresponds to
           // then call a function from a pointer table

           bondFunction* bonded = bondedInteractionFunctions[ftype]; 

           // compute all forces for ftype
           bonded(idef.iparams, idef.il[ftype], ...);
       }

       // reduce thread output
   }

|Gromacs| supports a lot of different listed interaction types, such as different
types of bonds, angles and proper and improper dihedrals. These different types
require different handling and finally the right force kernel chosen from a table
of function pointers.
The handling code required to correctly branch out to all the different cases
results in quite a deep call stack, a lot of branching logic and ends up accounting
for a fair part of the overall complexity, which should ideally just consist of
the type-specific force calculation implementations.


A type-aware approach to listed forces
--------------------------------------

NB-LIB aims to reduce the overall code complexity with a type-aware data format
where each interaction type is implemented as a separate (C++)-type.
The format for a given interaction type looks like this:

.. code:: cpp

   template <class Interaction>
   struct InteractionData
   {
       std::vector<Index<Interaction>> indices;
       std::vector<Interaction>        parameters;
   };

For each type of interaction, we store the interaction indices plus the
interaction parameters. While the (C++)-types are different, the actual data stored is
exactly the same: ``N+1`` integer indices per ``N``-center interaction plus the unique parameters.
An example for ``Interaction`` would be ``HarmonicBond``, the public part of which looks like this:

.. code:: cpp

   class HarmonicBond
   {
   public:
       // return lvalue ref for use with std::tie
       // in order to leverage std::tuple comparison ops
       const real& forceConstant();
       const real& equilDistance();
   };

The ``Index`` traits class deduces to ``std::array<int, 3>``, because
for each harmonic bond, we need two ``int``\ s for the coordinate
indices and a third ``int`` to look up the bond parameters in the
``parameters`` vector. For angles and dihedrals, the ``Index`` trait
would add an additional one or two ``int``\ s to hold the additional
coordinate indices.

Finally, we gather all types of interactions in a
``std::tuple``, such that the complete definition for listed forces
in NB-LIB looks like this:

.. code:: cpp

   using ListedInteractions = std::tuple<InteractionData<HarmonicBond>, ..., InteractionData<HarmonicAngle>, ...>;

One important property of ``ListedInteractions`` is that it stores exactly the same information as ``InteractionDefinitions``
and therefore conversion in either direction is easy to implement.


The NB-LIB listed forces pipeline
---------------------------------

Given the listed interaction data provided in the format described above,
the steps required to calculate the corresponding forces
are, in brief: 

  * Loop over all interaction types
  * Loop over all interactions for given type
  * Call interaction type kernel, store forces and return energy


This procedure is identical to the current implementation in |Gromacs|.
In actual code, the first step looks like this:

.. code:: cpp

   template<class Buffer, class Pbc>
   auto reduceListedForces(const ListedInteractions& interactions,
                           const std::vector<gmx::RVec>& x,
                           Buffer* forces,
                           const Pbc& pbc)
   {
       std::array<real, std::tuple_size<ListedInteractions>::value> energies;

       // lambda function, will be applied to each type
       auto computeForceType = [forces, &x, &energies, &pbc](const auto& ielem) {
           real energy = computeForces(ielem.indices, ielem.parameters, x, forces, pbc);
           energies[FindIndex<std::decay_t<decltype(ilem)>, ListedInteractions>{}] = energy;
       };

       // apply the lambda to all bond types
       for_each_tuple(computeForceType, interactions);

       return energies;
   }

With the help of a generic lambda and C++17’s ``std::apply`` in the
one-liner ``for_each_tuple``, we can generate the loop over the
different types in the tuple quite effortlessly. While
``reduceListedForces`` implements a loop over the interaction types, the
next layer, ``computeForces`` implements a loop over all interactions of
a given type:

.. code:: cpp

   template <class Index, class InteractionType, class Buffer, class Pbc>
   real computeForces(const std::vector<Index>& indices,
                      const std::vector<InteractionType>& iParams,
                      const std::vector<gmx::RVec>& x,
                      Buffer* forces,
                      const Pbc& pbc)
   {
       real Epot = 0.0;

       for (const auto& index : indices)
       {
           Epot += dispatchInteraction(index, iParams, x, forces);
       }

       return Epot;
   }

Compared to the union data type approach where this loop has been manually
implemented for all interaction types, in NB-LIB, only a single implementation
is required.

We’re now down to the level of individual bonds, angles and dihedrals.
At this point, the next steps depend on the actual type of the
interaction. But instead of dispatching each harmonic bond, cubic bond,
harmonic angle and so on to their seperate paths just yet, we just
differentiate based on the number of interaction centers for now.
Through overload resolution, the appropriate version
``dispatchInteraction`` gets called now, such as this one for the case
of 2-center interactions:

.. code:: cpp

   template <class Buffer, class TwoCenterType, class Pbc>
   std::enable_if_t<IsTwoCenter<TwoCenterType>::value, real>
   dispatchInteraction(const InteractionIndex<TwoCenterType>& index,
                       const std::vector<TwoCenterType>& bondInstances,
                       const std::vector<gmx::RVec>& x,
                       Buffer* forces,
                       const Pbc& pbc)
   {
       int i = std::get<0>(index);
       int j = std::get<1>(index);
       const gmx::RVec& x1 = x[i];
       const gmx::RVec& x2 = x[j];
       const TwoCenterType& bond = bondInstances[std::get<2>(index)];

       gmx::RVec dx;
       // calculate x1 - x2 modulo pbc
       pbc.dxAiuc(x1, x2, dx);
       real dr2 = dot(dx, dx);
       real dr  = std::sqrt(dr2);

       auto [force, energy] = bondKernel(dr, bond);

       // avoid division by 0
       if (dr2 != 0.0)
       {
           force /= dr;
           detail::spreadTwoCenterForces(force, dx, &(*forces)[i], &(*forces)[j]);
       }

       return energy;
   }

We can again observe that common parts among different 2-center interaction types
are reused. The common parts are 

 * coordinate retrieval
 * computation of the scalar distance
 * spreading of the scalar part of the force to the two centers

The only remaining thing to do now is to call the actual
kernel to compute the force. Since ``bond`` has a distinct type, we can
again use overload resolution:

.. code:: cpp

   template <class T>
   auto bondKernel(T dr, const HarmonicBond& bond)
   {
       return harmonicScalarForce(bond.forceConstant(), bond.equilDistance(), dr);
   }

and call the actual kernel, which in its simplest form for a harmonic
bond looks like this:

.. code:: cpp

   template <class T>
   std::tuple<T, T> harmonicScalarForce(T k, T x0, T x)
   {
       real dx  = x - x0;
       real dx2 = dx * dx;

       real force = -k * dx;
       real epot = 0.5 * k * dx2;

       return std::make_tuple(force, epot);

       /* That was 6 flops */
   }

That’s it! The approach outlined here manages to reuse (between different types)
a significant part of the code that feeds input data to force kernels.
Notably, not a single ``if(ftype)`` is required to implement the control flow.
The remaining parts for a feature complete implementation are
overloads of ``dispatchInteraction`` for the 3- to 5-center interactions and
the type-aware wrappers for all the different kernels implemented in
|Gromacs|. They have been omitted for brevity.

A note on **multithreading**: multithreading is handled above the top-level
``reduceListedForces`` described here. For parallel execution, the
input ``ListedInteractions`` tuple is split into ``nThreads`` parts and a
``Buffer`` object is set up for each thread. ``reduceListedForces`` is then
called once by each thread with the assigned fraction of ``ListedInteractions``
and the ``Buffer`` as argument.
The lifetime of the ``ListedInteractions`` splits is coupled to the domain decomposition.

Summary
-------

NB-LIB listed forces employs a (C++)-type aware data format that
is otherwise equivalent to its counter-part in |Gromacs|.
The type-aware data format is then used to simplify the "routing" layer that
connects data input to the appropriate kernels. Thanks to static branching and polymorphism,
increased code reuse and simplified branching logic could be achieved.
**The force kernels themselves do not need to be changed and NB-LIB refers to
|Gromacs| for their implementation.**


Outlook
-------

The data flow management for listed forces described here allows further
improvements to be implemented:

* Aggregate interaction types: fuse interactions of different types into
  aggregated types. For example, a dihedral interaction and the bonds and angles
  that are present among the same four particle indices can be combined into a single
  aggregated interaction. This allows to reuse the particle coordinates loaded from memory
  for multiple types and also combines the store operations for the forces.
  Type aggregates also likely simplify an efficient GPU implementation of listed forces.

* Separation of a topology containing both parameter sets for a system state A and B into two
  separate topologies for the A and B system states.
