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

#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/definitions.h"
#include "nblib/util/annotation.hpp"

namespace nblib
{

template<class I, class = void>
struct HasCharge : util::integral_constant<int, 0>
{
};

template<class I>
struct HasCharge<I, std::enable_if_t<Contains<I, ChargedListedTypes>{}>> :
    util::integral_constant<int, 1>
{
};

template<class I, class = void>
struct HasTwoCenterAggregate : util::integral_constant<int, 0>
{
};

template<class I>
struct HasTwoCenterAggregate<I, std::void_t<typename I::TwoCenterAggregateType>> :
    util::integral_constant<int, 1>
{
};

template<class I, class = void>
struct HasThreeCenterAggregate : util::integral_constant<int, 0>
{
};

template<class I>
struct HasThreeCenterAggregate<I, std::void_t<typename I::ThreeCenterAggregateType>> :
    util::integral_constant<int, 1>
{
};

template<class I, class = void>
struct HasPairAggregate : util::integral_constant<int, 0>
{
};

template<class I>
struct HasPairAggregate<I, std::void_t<typename I::PairAggregateType>> : util::integral_constant<int, 1>
{
};

template<class I, class = void>
struct IsAggregate : util::integral_constant<int, 0>
{
};

template<class I>
struct IsAggregate<I, std::void_t<typename I::CarrierType>> : util::integral_constant<int, 1>
{
};

//! \internal \brief determines the energy storage location of the carrier part for InteractionTypes without aggregates
template<class InteractionType, class = void>
struct CarrierIndex : util::integral_constant<size_t, FindIndex<InteractionType, AllListedTypes>{}>
{
};

//! \internal \brief determines the energy storage location of the carrier part for InteractionTypes with aggregates
template<class InteractionType>
struct CarrierIndex<InteractionType, std::void_t<typename InteractionType::CarrierType>> :
    util::integral_constant<size_t, FindIndex<typename InteractionType::CarrierType, AllListedTypes>{}>
{
};

//! \internal \brief determines the energy storage location of the 2-C aggregate part for InteractionTypes without aggregates
template<class InteractionType, class = void>
struct TwoCenterAggregateIndex : util::integral_constant<size_t, 0>
{
};

//! \internal \brief determines the energy storage location of the 2-C aggregate part for InteractionTypes with 2-C aggregates
template<class InteractionType>
struct TwoCenterAggregateIndex<InteractionType, std::void_t<typename InteractionType::TwoCenterAggregateType>> :
    util::integral_constant<size_t, FindIndex<typename InteractionType::TwoCenterAggregateType, AllListedTypes>{}>
{
};

//! \internal \brief determines the energy storage location of the 3-C aggregate part for InteractionTypes without aggregates
template<class InteractionType, class = void>
struct ThreeCenterAggregateIndex : util::integral_constant<size_t, 0>
{
};

//! \internal \brief determines the energy storage location of the 3-C aggregate part for InteractionTypes with 3-C aggregates
template<class InteractionType>
struct ThreeCenterAggregateIndex<InteractionType, std::void_t<typename InteractionType::ThreeCenterAggregateType>> :
    util::integral_constant<size_t, FindIndex<typename InteractionType::ThreeCenterAggregateType, AllListedTypes>{}>
{
};

/*! \brief array used to store potential energies of all listed interaction types
 *
 * In addition to all the listed interaction types, the last 3 elements are eVdW, eCoulomb, and dvdl
 */
using ListedEnergies = util::array<real, TypeListSize<AllListedTypes>{} + 3>;

struct VdwIndex : util::integral_constant<int, TypeListSize<AllListedTypes>{}>
{
};

struct CoulombIndex : util::integral_constant<int, TypeListSize<AllListedTypes>{} + 1>
{
};

struct FepIndex : util::integral_constant<int, TypeListSize<AllListedTypes>{} + 2>
{
};

struct GpuListedEnergySize : util::integral_constant<int, TypeListSize<GpuListedTypes>{} + 3>
{
};

struct GpuVdwIndex : util::integral_constant<int, TypeListSize<GpuListedTypes>{}>
{
};

struct GpuCoulombIndex : util::integral_constant<int, TypeListSize<GpuListedTypes>{} + 1>
{
};

struct GpuFepIndex : util::integral_constant<int, TypeListSize<GpuListedTypes>{} + 2>
{
};

} // namespace nblib
#endif // NBLIB_LISTEDFORCES_TRAITS_H
