/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * This implements conversion utilities between the internal
 * representations of the listed forces parameters for NBLIB
 * and that of the GROMACS backend
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTEDFORCES_CONVERSION_HPP
#define NBLIB_LISTEDFORCES_CONVERSION_HPP

#include <memory>

#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "nblib/listed_forces/traits.h"

namespace nblib
{

/*! this trait maps nblib InteractionTypes to the corresponding gmx enum
 *
 * \tparam InteractionType
 */
template<class InteractionType>
struct ListedIndex
{
};

template<>
struct ListedIndex<HarmonicBondType> : std::integral_constant<int, F_BONDS>
{
};

template<>
struct ListedIndex<MorseBondType> : std::integral_constant<int, F_MORSE>
{
};

template<>
struct ListedIndex<FENEBondType> : std::integral_constant<int, F_FENEBONDS>
{
};

template<>
struct ListedIndex<CubicBondType> : std::integral_constant<int, F_CUBICBONDS>
{
};

template<>
struct ListedIndex<G96BondType> : std::integral_constant<int, F_G96BONDS>
{
};

template<>
struct ListedIndex<HarmonicAngle> : std::integral_constant<int, F_ANGLES>
{
};

template<>
struct ListedIndex<ProperDihedral> : std::integral_constant<int, F_PDIHS>
{
};

template<>
struct ListedIndex<ImproperDihedral> : std::integral_constant<int, F_IDIHS>
{
};

template<>
struct ListedIndex<RyckaertBellemanDihedral> : std::integral_constant<int, F_RBDIHS>
{
};

template<class InteractionType>
struct ListedTypeIsImplemented : std::bool_constant<detail::HasValueMember<ListedIndex<InteractionType>>::value>
{
};

namespace detail
{

template<class InteractionData>
void transferParameters([[maybe_unused]] const InteractionData& interactionData,
                        [[maybe_unused]] gmx_ffparams_t&        gmx_params)
{
}

//! \brief Harmonic bond parameter conversion function
template<>
void transferParameters(const ListedTypeData<HarmonicBondType>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& hbond : interactions.parameters)
    {
        t_iparams param;
        param.harmonic.krA = hbond.forceConstant();
        param.harmonic.rA  = hbond.equilConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Harmonic angle parameter conversion function
template<>
void transferParameters(const ListedTypeData<HarmonicAngle>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& angle : interactions.parameters)
    {
        t_iparams param;
        param.harmonic.krA = angle.forceConstant();
        param.harmonic.rA  = angle.equilConstant() / DEG2RAD;
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Proper dihedral parameter conversion function
template<>
void transferParameters(const ListedTypeData<ProperDihedral>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& dihedral : interactions.parameters)
    {
        t_iparams param;
        param.pdihs.phiA = dihedral.equilDistance() / DEG2RAD;
        param.pdihs.cpA  = dihedral.forceConstant();
        param.pdihs.mult = dihedral.multiplicity();
        gmx_params.iparams.push_back(param);
    }
}

template<class TwoCenterType>
std::enable_if_t<Contains<TwoCenterType, SupportedTwoCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<TwoCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int parameterIndex = index[2] + offset;
        idef.il[ListedIndex<TwoCenterType>::value].iatoms.push_back(parameterIndex);
        idef.il[ListedIndex<TwoCenterType>::value].iatoms.push_back(index[0]);
        idef.il[ListedIndex<TwoCenterType>::value].iatoms.push_back(index[1]);
    }
}

template<class ThreeCenterType>
std::enable_if_t<Contains<ThreeCenterType, SupportedThreeCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<ThreeCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int parameterIndex = index[3] + offset;
        idef.il[ListedIndex<ThreeCenterType>::value].iatoms.push_back(parameterIndex);
        idef.il[ListedIndex<ThreeCenterType>::value].iatoms.push_back(index[0]);
        idef.il[ListedIndex<ThreeCenterType>::value].iatoms.push_back(index[1]);
        idef.il[ListedIndex<ThreeCenterType>::value].iatoms.push_back(index[2]);
    }
}

template<class FourCenterType>
std::enable_if_t<Contains<FourCenterType, SupportedFourCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<FourCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int parameterIndex = index[4] + offset;
        idef.il[ListedIndex<FourCenterType>::value].iatoms.push_back(parameterIndex);
        idef.il[ListedIndex<FourCenterType>::value].iatoms.push_back(index[0]);
        idef.il[ListedIndex<FourCenterType>::value].iatoms.push_back(index[1]);
        idef.il[ListedIndex<FourCenterType>::value].iatoms.push_back(index[2]);
        idef.il[ListedIndex<FourCenterType>::value].iatoms.push_back(index[3]);
    }
}

template<class FiveCenterType>
std::enable_if_t<Contains<FiveCenterType, SupportedFiveCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<FiveCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int parameterIndex = index[5] + offset;
        idef.il[ListedIndex<FiveCenterType>::value].iatoms.push_back(parameterIndex);
        idef.il[ListedIndex<FiveCenterType>::value].iatoms.push_back(index[0]);
        idef.il[ListedIndex<FiveCenterType>::value].iatoms.push_back(index[1]);
        idef.il[ListedIndex<FiveCenterType>::value].iatoms.push_back(index[2]);
        idef.il[ListedIndex<FiveCenterType>::value].iatoms.push_back(index[3]);
        idef.il[ListedIndex<FiveCenterType>::value].iatoms.push_back(index[4]);
    }
}

template<template<class> class Container, class InteractionType>
void transferIndices(const Container<InteractionType>&  interactionData,
                     InteractionDefinitions& idef,
                     [[maybe_unused]] int offset)
{
    if constexpr (ListedTypeIsImplemented<InteractionType>{})
    {
        transferIndicesImpl(interactionData, idef, offset);
    }
}

} // namespace detail

/*! \brief function to convert from nblib-ListedInteractionData to gmx-InteractionDefinitions
 *
 * Currently only supports harmonic bonds and angles, other types are ignored
 *
 * \param interactions
 * \return
 */
std::tuple<std::unique_ptr<InteractionDefinitions>, std::unique_ptr<gmx_ffparams_t>>
createFFparams(const ListedInteractionData& interactions);

std::tuple<std::unique_ptr<InteractionDefinitions>, std::unique_ptr<gmx_ffparams_t>>
createFFparams(const ListedInteractionData& interactions)
{
    std::unique_ptr<gmx_ffparams_t> ffparamsHolder = std::make_unique<gmx_ffparams_t>();
    std::unique_ptr<InteractionDefinitions> idefHolder = std::make_unique<InteractionDefinitions>(*ffparamsHolder);

    gmx_ffparams_t& ffparams = *ffparamsHolder;
    InteractionDefinitions& idef = *idefHolder;

    auto copyParamsOneType = [&ffparams](const auto& interactionElement)
    {
        detail::transferParameters(interactionElement, ffparams);
    };
    for_each_tuple(copyParamsOneType, interactions);

    // since gmx_ffparams_t.iparams is a flattened vector over all interaction types,
    // we need to compute offsets for each type to know where the parameters for each type start
    // in the flattened iparams vectors
    int acc = 0;
    std::array<int, std::tuple_size_v<ListedInteractionData>> indexOffsets{0};
    auto extractNIndices = [&indexOffsets, &acc](const auto& interactionElement)
    {
        constexpr int elementIndex = FindIndex<std::decay_t<decltype(interactionElement)>, ListedInteractionData>::value;
        indexOffsets[elementIndex] = acc;
        acc += interactionElement.parameters.size();
    };
    for_each_tuple(extractNIndices, interactions);

    auto copyIndicesOneType = [&idef, &indexOffsets](const auto& interactionElement)
    {
        constexpr int elementIndex = FindIndex<std::decay_t<decltype(interactionElement)>, ListedInteractionData>::value;
        detail::transferIndices(interactionElement, idef, indexOffsets[elementIndex]);
    };
    for_each_tuple(copyIndicesOneType, interactions);

    return std::make_tuple(std::move(idefHolder), std::move(ffparamsHolder));
}


} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CONVERSION_HPP
