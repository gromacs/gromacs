Adding New Listed-Interaction Types in NB-LIB
=============================================

NB-LIB currently has code paths for listed interactions that occur between two, three, four and five different particles.
It is easy to extend NB-LIB to support novel formulations of particle interactions by modifying the following three files.

Two center interactions must use the distance between the centers as an input to the force kernel.
Three center interactions take the form ``(particleI, particleJ, ParticleK)``.
In this case, the middle particle, ``particleJ`` is taken as the center around which the angle is computed.
This angle must be an input to a three center force kernel.
Likewise for four center interactions, the dihedral angle phi must be an input to the force kernel.
Accepting these constraints, it is possible to add a new kernel by modifying the following three files.

1) bondtypes.h_
2) definitions.h_
3) kernels.hpp_

.. _bondtypes.h:

1) bondtypes.h
---------------

This file contains one ``struct`` for each interaction type parameter set.
New interaction types are added here as separate structs. There
are no content requirements, but for convenience, the existing ``NAMED_MEBERS``
macro in combination with inheriting from a ``std::tuple`` or ``std::array``
may be used. The macro can be used to define the
parameter names for the corresponding setters and getters.
For example, ``NAMED_MEMBERS(forceConstant, equilDistance)`` will expand to

.. code:: cpp

   inline auto& forceConstant() { return std::get<0>(*this); }
   inline auto& equilDistance() { return std::get<1>(*this); }
   inline const auto& forceConstant() const { return std::get<0>(*this); }
   inline const auto& equilDistance() const { return std::get<1>(*this); }

Putting everything together, one could define the complete parameter set for a new interaction type as follows.

.. code:: cpp

   /*! \brief new bond type
    *
    * V(r; forceConstant, equilDistance, scaleFactor)
    *      = forceConstant * exp( (r - equilDistance) / scaleFactor)
    */
   struct NewBondType : public std::tuple<real, real, int>
   {
       NewBondType() = default;
       NewBondType(ForceConstant f, EquilDistance d, ScaleFactor s) :
           std::tuple<real, real, int>{ f, d, s }
       {
       }

       NAMED_MEMBERS(forceConstant, equilDistance, scaleFactor)
   };

.. _definitions.h:

2) definitions.h
------------------------

This file begins with pre-processor macro lists that classify concrete interaction types into two, three, four and five center types.
To add a new type, the user must add the interaction type parameter struct name to the macro of the correct center number.
In this case, ``NewBondType`` is an example of a two center interaction.
As such it would get added to the ``SUPPORTED_TWO_CENTER_TYPES`` macro.
Assuming that the only other two center interaction is called ``DefaultBond``, the result would look like the following snippet.

.. code:: cpp

    #define SUPPORTED_TWO_CENTER_TYPES DefaultBond, NewBondType

.. _kernels.hpp:

Adding ``NewBondType`` to this macro ensures that the NB-LIB ``molecule``
class ``addInteraction`` function supports adding the new bond type
and includes it in the listed interaction data that the ``topology`` class
provides.

Note that, as of C++17, there's no alternative to preprocessor macros for adding
the required template instantiations controlled through the macros described here.
In NB-LIB, the design decision we took, was that we did not want to expose a templated
interface in a user header and it is for this reason that we explicitly need
to instantiate the interface with all the supported listed interaction types defined
in this macro.

3) kernels.hpp
---------------------

In this file the actual force kernels for each interaction type are implemented.
Each kernel call is templated to allow various precisions and is
accessed through an overload ``bondKernel`` that extracts the relevant
parameters from a ``const NewBondType&`` argument.
The kernel return type is always an ``std::tuple`` of the force and the potential.

.. code:: cpp

   /*! \brief kernel to calculate the new bond type force
    *
    * \param k     Force constant
    * \param x0    Equilibrium distance
    * \param scale The scaling factor
    * \param x     Input bond length
    *
    * \return tuple<force, potential energy>
    */
   template <class T>
   std::tuple<T, T> newBondForce(T k, T x0, T scale, T x)
   {
       real exponent = std::exp( (x - x0) / scale);
       real epot = k * exponent;
       real force =  epot / scale;
       return std::make_tuple(force, epot);
   }

  template <class T>
  inline auto bondKernel(T dr, const NewBondType& bond)
  {
      return newBondForce(bond.forceConstant(), bond.equilDistance(), bond.scaleFactor(), dr);
  }

