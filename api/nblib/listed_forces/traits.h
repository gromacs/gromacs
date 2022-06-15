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
 * These traits defined here for supported nblib listed interaction data types
 * are used to control the dataflow in dataflow.hpp
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_LISTEDFORCES_TRAITS_H
#define NBLIB_LISTEDFORCES_TRAITS_H

#include <numeric>

#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/definitions.h"

namespace nblib
{

template<class I, class = void>
struct HasTwoCenterAggregate : std::false_type
{
};

template<class I>
struct HasTwoCenterAggregate<I, std::void_t<typename I::TwoCenterAggregateType>> : std::true_type
{
};

template<class I, class = void>
struct HasThreeCenterAggregate : std::false_type
{
};

template<class I>
struct HasThreeCenterAggregate<I, std::void_t<typename I::ThreeCenterAggregateType>> : std::true_type
{
};

//! \internal \brief determines the energy storage location of the carrier part for InteractionTypes without aggregates
template<class InteractionType, class = void>
struct CarrierIndex :
    std::integral_constant<size_t, FindIndex<InteractionType, ListedInteractionData>{}>
{
};

//! \internal \brief determines the energy storage location of the carrier part for InteractionTypes with aggregates
template<class InteractionType>
struct CarrierIndex<InteractionType, std::void_t<typename InteractionType::CarrierType>> :
    std::integral_constant<size_t, FindIndex<typename InteractionType::CarrierType, ListedInteractionData>{}>
{
};

//! \internal \brief determines the energy storage location of the 2-C aggregate part for InteractionTypes without aggregates
template<class InteractionType, class = void>
struct TwoCenterAggregateIndex : std::integral_constant<size_t, 0>
{
};

//! \internal \brief determines the energy storage location of the 2-C aggregate part for InteractionTypes with 2-C aggregates
template<class InteractionType>
struct TwoCenterAggregateIndex<InteractionType, std::void_t<typename InteractionType::TwoCenterAggregateType>> :
    std::integral_constant<size_t, FindIndex<typename InteractionType::TwoCenterAggregateType, ListedInteractionData>{}>
{
};

//! \internal \brief determines the energy storage location of the 3-C aggregate part for InteractionTypes without aggregates
template<class InteractionType, class = void>
struct ThreeCenterAggregateIndex : std::integral_constant<size_t, 0>
{
};

//! \internal \brief determines the energy storage location of the 3-C aggregate part for InteractionTypes with 3-C aggregates
template<class InteractionType>
struct ThreeCenterAggregateIndex<InteractionType, std::void_t<typename InteractionType::ThreeCenterAggregateType>> :
    std::integral_constant<size_t, FindIndex<typename InteractionType::ThreeCenterAggregateType, ListedInteractionData>{}>
{
};

/*! \brief return type to hold the energies of the different overloads of "dispatchInteraction"
 * \internal
 *
 * \tparam T
 */
template<class T>
class KernelEnergy
{
public:
    KernelEnergy() : energies_{ 0, 0, 0, 0 } {}

    T&       carrier() { return energies_[0]; }
    const T& carrier() const { return energies_[0]; }

    T&       twoCenterAggregate() { return energies_[1]; }
    const T& twoCenterAggregate() const { return energies_[1]; }

    T&       threeCenterAggregate() { return energies_[2]; }
    const T& threeCenterAggregate() const { return energies_[2]; }

    T&       freeEnergyDerivative() { return energies_[3]; }
    const T& freeEnergyDerivative() const { return energies_[3]; }

    KernelEnergy& operator+=(const KernelEnergy& other)
    {
        for (size_t i = 0; i < energies_.size(); ++i)
        {
            energies_[i] += other.energies_[i];
        }
        return *this;
    }

    operator T() const { return std::accumulate(begin(energies_), end(energies_), T{}); }

private:
    std::array<T, 4> energies_;
};

template<class BasicVector>
using BasicVectorValueType_t = std::remove_all_extents_t<typename BasicVector::RawArray>;

} // namespace nblib
#endif // NBLIB_LISTEDFORCES_TRAITS_H
