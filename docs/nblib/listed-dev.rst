Adding New Listed-Interaction Types in NB-LIB
=============================================

NB-LIB currently has code paths for listed interactions that occur between two, three, four and five different particles.
To extend NB-LIB to support more types of particle interactions, modify the following three files.

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

This file contains one C++ type to store the parameters for each interaction type.
New interaction types are added here as separate C++ types.
The interface of these types is completely unrestricted.
The only requirements are equality and less than comparison, and that the interface be
compatible with the corresponding (user-added) kernel.

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

Adding ``NewBondType`` to this macro ensures that the NBLIB ``molecule``
class ``addInteraction`` function supports adding the new bond type
and includes it in the listed interaction data that the ``topology`` class
provides. The ``SUPPORTED_TWO_CENTER_TYPES`` macro is immediately converted into a
C++ type list that is implemented as a variadic template. The type list
is then used to define all the dependent data structures. Apart from creating
the type list, the only place where the macro is needed is explicit template instantiation.

Note that, as of C++17, there's no alternative to preprocessor macros for adding
the required template instantiations controlled through the macros described here.
(Other than manually adding the template instantiations, which would require the instantiation list
of several templates to be updated each time a new interaction type is added. Compared to the preprocessor
based solution where just a single macro has to be extended, this would clearly be an inferior solution.)
In NBLIB, the design decision we took, was that we did not want to expose a templated
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
  inline std::tuple<T, T> bondKernel(T dr, const NewBondType& bond)
  {
      return newBondForce(bond.forceConstant(), bond.equilDistance(), bond.scaleFactor(), dr);
  }

