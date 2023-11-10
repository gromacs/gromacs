/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements load/store operations for host and device as well as
 * selectors tags to differentiate between FEP yes/no shift forces yes/no.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTEDFORCES_SELECTORS_HPP
#define NBLIB_LISTEDFORCES_SELECTORS_HPP

#include <type_traits>
#include <utility>

#include "nblib/listed_forces/definitions.h"
#include "nblib/util/array.hpp"

namespace nblib
{

/*! \brief 3D vector type for operations on local stack variables
 *
 * The 3D-Vector type used to call the dispatchInteraction functions depend on the
 * platform and may be 3- or 4-elements wide and may be StrongType-tagged for overload resolution,
 * e.g. atomic stores for GPU force buffers. Instead of argument type deduction we use this
 * as a fixed type for local stack variables to keep load/store overload resolution clean and simple.
 */
template<class T>
using StackVec3 = util::array<T, 3>;

//! \brief tag to mark non-FEP code-paths
struct NoFepLambdaType
{
};

//! \brief a tag used to distinguish host from device-resident memory
template<class T>
using DeviceTag = StrongType<T, struct DeviceTagParameter>;

template<class ShiftForce, class = void>
struct HaveShifts : util::integral_constant<int, 1> // integral_constant for device compatibility
{
};

template<class ShiftForce>
struct HaveShifts<ShiftForce, std::enable_if_t<std::is_same_v<std::decay_t<ShiftForce>, std::nullptr_t>>> :
    util::integral_constant<int, 0>
{
};

//! \brief load B parameters if a real-valued lambda is passed
template<class Lambda, class InteractionType, std::enable_if_t<!std::is_same_v<Lambda, NoFepLambdaType>, int> = 0>
HOST_DEVICE_FUN inline InteractionType loadInteractionParameters(const InteractionType* parameters, int index)
{
    return parameters[index];
}

//! \brief send back empty parameters when not calculating FEP (will be discarded by the compiler)
template<class Lambda, class InteractionType, std::enable_if_t<std::is_same_v<Lambda, NoFepLambdaType>, int> = 0>
HOST_DEVICE_FUN inline InteractionType loadInteractionParameters(const InteractionType* /*parameters*/,
                                                                 int /*index*/)
{
    return InteractionType{};
}

//! \brief load a 3D XYZ vector from memory, input may be 4-wide
template<class MemVector>
HOST_DEVICE_FUN inline StackVec3<VectorValueType_t<MemVector>> loadVec(MemVector vector)
{
    using T = VectorValueType_t<MemVector>;
    return StackVec3<T>{ vector[0], vector[1], vector[2] };
}

//! \brief store forces back on the host
template<class MemVector, class T>
inline void addForce(MemVector* location, const StackVec3<T>& force)
{
    (*location)[0] += force[0];
    (*location)[1] += force[1];
    (*location)[2] += force[2];
}

//! \brief store forces back to global device memory
template<class MemVector, class T>
DEVICE_FUN inline void addForce(DeviceTag<MemVector>* location, const StackVec3<T>& force)
{
    atomicAdd(&location->value()[0], force[0]);
    atomicAdd(&location->value()[1], force[1]);
    atomicAdd(&location->value()[2], force[2]);
}

//! \brief no-op overload, used in the absence of shift forces
template<class T>
inline void addForce(std::nullptr_t*, const StackVec3<T>&)
{
}

//! \brief no-op overload, used in the absence of shift forces
template<class T>
DEVICE_FUN inline void addForce(DeviceTag<std::nullptr_t>*, const StackVec3<T>&)
{
}

/*! \brief return type to hold the energies of the different overloads of "dispatchInteraction"
 * \internal
 *
 * \tparam T
 */
template<class T>
class KernelEnergy
{
public:
    HOST_DEVICE_FUN KernelEnergy() : energies_{ 0, 0, 0, 0, 0, 0 } {}

    HOST_DEVICE_FUN T&    carrier() { return energies_[0]; }
    HOST_DEVICE_FUN const T& carrier() const { return energies_[0]; }

    HOST_DEVICE_FUN T&    twoCenterAggregate() { return energies_[1]; }
    HOST_DEVICE_FUN const T& twoCenterAggregate() const { return energies_[1]; }

    HOST_DEVICE_FUN T&    threeCenterAggregate() { return energies_[2]; }
    HOST_DEVICE_FUN const T& threeCenterAggregate() const { return energies_[2]; }

    HOST_DEVICE_FUN T&    freeEnergyDerivative() { return energies_[3]; }
    HOST_DEVICE_FUN const T& freeEnergyDerivative() const { return energies_[3]; }

    HOST_DEVICE_FUN T&    eVdw() { return energies_[4]; }
    HOST_DEVICE_FUN const T& eVdw() const { return energies_[4]; }

    HOST_DEVICE_FUN T&    eCoul() { return energies_[5]; }
    HOST_DEVICE_FUN const T& eCoul() const { return energies_[5]; }

    HOST_DEVICE_FUN KernelEnergy& operator+=(const KernelEnergy& other)
    {
        energies_ += other.energies_;

        return *this;
    }

    HOST_DEVICE_FUN T potentialEnergy() const
    {
        return carrier() + twoCenterAggregate() + threeCenterAggregate();
    }

private:
    util::array<T, 6> energies_;
};

template<class InteractionType, class T, class Tuple, std::enable_if_t<Contains<InteractionType, PolarizationTypes>{}, int> = 0>
HOST_DEVICE_FUN inline void addEnergy(KernelEnergy<T>* energy, const Tuple& energyTuple)
{
    energy->carrier() += util::get<0>(energyTuple);
    energy->freeEnergyDerivative() += util::get<1>(energyTuple);
}

template<class InteractionType, class T, class Tuple, std::enable_if_t<Contains<InteractionType, PairListedTypes>{}, int> = 0>
HOST_DEVICE_FUN inline void addEnergy(KernelEnergy<T>* energy, const Tuple& energyTuple)
{
    energy->eVdw() += util::get<0>(energyTuple);
    energy->eCoul() += util::get<1>(energyTuple);
    energy->freeEnergyDerivative() += util::get<2>(energyTuple);
}

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_SELECTORS_HPP
