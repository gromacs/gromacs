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

#include <cmath>
#include <cstddef>

#include <array>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "listed_forces/conversionscommon.h"

#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"

#include "nblib/basicdefinitions.h"
#include "nblib/exception.h"
#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/definitions.h"
#include "nblib/util/traits.hpp"
#include "nblib/util/util.hpp"

namespace nblib
{

namespace detail
{

template<class InteractionData>
inline void transferParameters(const InteractionData& /* interactionDataA */,
                               const InteractionData& /* interactionDataB */,
                               gmx_ffparams_t& /* gmx_params */)
{
}

template<class InteractionData>
inline void transferParameters([[maybe_unused]] const InteractionData& interactionData,
                               [[maybe_unused]] gmx_ffparams_t&        gmx_params)
{
    transferParameters(interactionData, interactionData, gmx_params);
}

//! \brief Harmonic bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<HarmonicBondType>& interactionsA,
                               const ListedTypeData<HarmonicBondType>& interactionsB,
                               gmx_ffparams_t&                         gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("Harmonic bond interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactionsA.parameters[i].forceConstant();
        param.harmonic.rA  = interactionsA.parameters[i].equilConstant();
        param.harmonic.krB = interactionsB.parameters[i].forceConstant();
        param.harmonic.rB  = interactionsB.parameters[i].equilConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief G96 bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<G96BondType>& interactionsA,
                               const ListedTypeData<G96BondType>& interactionsB,
                               gmx_ffparams_t&                    gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("G96 bond interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactionsA.parameters[i].forceConstant();
        param.harmonic.rA  = interactionsA.parameters[i].equilConstant();
        param.harmonic.krB = interactionsB.parameters[i].forceConstant();
        param.harmonic.rB  = interactionsB.parameters[i].equilConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief FENE bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<FENEBondType>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& bond : interactions.parameters)
    {
        t_iparams param;
        param.fene.kb = bond.forceConstant();
        param.fene.bm = bond.equilConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Cubic bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<CubicBondType>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& bond : interactions.parameters)
    {
        t_iparams param;
        param.cubic.kcub = bond.cubicForceConstant();
        param.cubic.kb   = bond.quadraticForceConstant();
        param.cubic.b0   = bond.equilDistance();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Morse bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<MorseBondType>& interactionsA,
                               const ListedTypeData<MorseBondType>& interactionsB,
                               gmx_ffparams_t&                      gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("Morse bond interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.morse.b0A   = interactionsA.parameters[i].equilDistance();
        param.morse.betaA = interactionsA.parameters[i].exponent();
        param.morse.cbA   = interactionsA.parameters[i].forceConstant();
        param.morse.b0B   = interactionsB.parameters[i].equilDistance();
        param.morse.betaB = interactionsB.parameters[i].exponent();
        param.morse.cbB   = interactionsB.parameters[i].forceConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Pair LJ 1-4 interaction parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<PairLJType>& interactionsA,
                               const ListedTypeData<PairLJType>& interactionsB,
                               gmx_ffparams_t&                   gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("LJ1-4 pair interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.lj14.c6A  = interactionsA.parameters[i].c6();
        param.lj14.c12A = interactionsA.parameters[i].c12();
        param.lj14.c6B  = interactionsB.parameters[i].c6();
        param.lj14.c12B = interactionsB.parameters[i].c12();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Harmonic angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<HarmonicAngle>& interactionsA,
                               const ListedTypeData<HarmonicAngle>& interactionsB,
                               gmx_ffparams_t&                      gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("Harmonic angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactionsA.parameters[i].forceConstant();
        param.harmonic.rA  = interactionsA.parameters[i].equilConstant() / DEG2RAD;
        param.harmonic.krB = interactionsB.parameters[i].forceConstant();
        param.harmonic.rB  = interactionsB.parameters[i].equilConstant() / DEG2RAD;
        gmx_params.iparams.push_back(param);
    }
}

//! \brief G96 angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<G96Angle>& interactionsA,
                               const ListedTypeData<G96Angle>& interactionsB,
                               gmx_ffparams_t&                 gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("G96 angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactionsA.parameters[i].forceConstant();
        param.harmonic.rA  = interactionsA.parameters[i].equilConstant();
        param.harmonic.krB = interactionsB.parameters[i].forceConstant();
        param.harmonic.rB  = interactionsB.parameters[i].equilConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Linear angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<LinearAngle>& interactionsA,
                               const ListedTypeData<LinearAngle>& interactionsB,
                               gmx_ffparams_t&                    gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("Linear angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.linangle.klinA = interactionsA.parameters[i].forceConstant();
        param.linangle.aA    = interactionsA.parameters[i].equilConstant();
        param.linangle.klinB = interactionsB.parameters[i].forceConstant();
        param.linangle.aB    = interactionsB.parameters[i].equilConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Restricted angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<RestrictedAngle>& interactionsA,
                               const ListedTypeData<RestrictedAngle>& interactionsB,
                               gmx_ffparams_t&                        gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("Restricted angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactionsA.parameters[i].forceConstant();
        param.harmonic.rA  = (std::acos(interactionsA.parameters[i].equilConstant())) / DEG2RAD;
        param.harmonic.krB = interactionsB.parameters[i].forceConstant();
        param.harmonic.rB  = (std::acos(interactionsB.parameters[i].equilConstant())) / DEG2RAD;
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Quartic angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<QuarticAngle>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& angle : interactions.parameters)
    {
        t_iparams param;
        param.qangle.theta = angle.equilConstant() / DEG2RAD;
        param.qangle.c[0]  = angle.forceConstant(0);
        param.qangle.c[1]  = angle.forceConstant(1);
        param.qangle.c[2]  = angle.forceConstant(2);
        param.qangle.c[3]  = angle.forceConstant(3);
        param.qangle.c[4]  = angle.forceConstant(4);
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Cross Bond-Bond interaction parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<CrossBondBond>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& bond : interactions.parameters)
    {
        t_iparams param;
        param.cross_bb.krr = bond.forceConstant();
        param.cross_bb.r1e = bond.equilDistanceIJ();
        param.cross_bb.r2e = bond.equilDistanceKJ();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Cross Bond-Angle interaction parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<CrossBondAngle>& interactions, gmx_ffparams_t& gmx_params)
{
    for (const auto& bond : interactions.parameters)
    {
        t_iparams param;
        param.cross_ba.krt = bond.forceConstant();
        param.cross_ba.r1e = bond.equilDistanceIJ();
        param.cross_ba.r2e = bond.equilDistanceKJ();
        param.cross_ba.r3e = bond.equilDistanceIK();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Proper dihedral parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<ProperDihedral>& interactionsA,
                               const ListedTypeData<ProperDihedral>& interactionsB,
                               gmx_ffparams_t&                       gmx_params)
{
    if (interactionsA.parameters.size() != interactionsB.parameters.size()
        && interactionsA.indices.size() != interactionsB.indices.size())
    {
        throw InputException("Proper dihedral interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.pdihs.phiA = interactionsA.parameters[i].equilDistance() / DEG2RAD;
        param.pdihs.cpA  = interactionsA.parameters[i].forceConstant();
        param.pdihs.mult = interactionsA.parameters[i].multiplicity();
        param.pdihs.phiB = interactionsB.parameters[i].equilDistance() / DEG2RAD;
        param.pdihs.cpB  = interactionsB.parameters[i].forceConstant();
        gmx_params.iparams.push_back(param);
    }
}

//! \brief Ryckaert-Belleman dihedral parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<RyckaertBellemanDihedral>& interactionsA,
                               const ListedTypeData<RyckaertBellemanDihedral>& interactionsB,
                               gmx_ffparams_t&                                 gmx_params)
{
    for (size_t i = 0; i < interactionsA.parameters.size(); i++)
    {
        t_iparams param;
        param.rbdihs.rbcA[0] = interactionsA.parameters[i][0];
        param.rbdihs.rbcA[1] = interactionsA.parameters[i][1];
        param.rbdihs.rbcA[2] = interactionsA.parameters[i][2];
        param.rbdihs.rbcA[3] = interactionsA.parameters[i][3];
        param.rbdihs.rbcA[4] = interactionsA.parameters[i][4];
        param.rbdihs.rbcA[5] = interactionsA.parameters[i][5];
        param.rbdihs.rbcB[0] = interactionsB.parameters[i][0];
        param.rbdihs.rbcB[1] = interactionsB.parameters[i][1];
        param.rbdihs.rbcB[2] = interactionsB.parameters[i][2];
        param.rbdihs.rbcB[3] = interactionsB.parameters[i][3];
        param.rbdihs.rbcB[4] = interactionsB.parameters[i][4];
        param.rbdihs.rbcB[5] = interactionsB.parameters[i][5];
        gmx_params.iparams.push_back(param);
    }
}

template<class TwoCenterType>
inline std::enable_if_t<Contains<TwoCenterType, SupportedTwoCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<TwoCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int  parameterIndex = index[2] + offset;
        auto gmxListedID    = FindIndex<TwoCenterType, GmxToNblibMapping>::value;
        idef.il[gmxListedID].iatoms.push_back(parameterIndex);
        idef.il[gmxListedID].iatoms.push_back(index[0]);
        idef.il[gmxListedID].iatoms.push_back(index[1]);
    }
}

template<class ThreeCenterType>
inline std::enable_if_t<Contains<ThreeCenterType, SupportedThreeCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<ThreeCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int  parameterIndex = index[3] + offset;
        auto gmxListedID    = FindIndex<ThreeCenterType, GmxToNblibMapping>::value;
        idef.il[gmxListedID].iatoms.push_back(parameterIndex);
        idef.il[gmxListedID].iatoms.push_back(index[0]);
        idef.il[gmxListedID].iatoms.push_back(index[1]);
        idef.il[gmxListedID].iatoms.push_back(index[2]);
    }
}

template<class FourCenterType>
inline std::enable_if_t<Contains<FourCenterType, SupportedFourCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<FourCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int  parameterIndex = index[4] + offset;
        auto gmxListedID    = FindIndex<FourCenterType, GmxToNblibMapping>::value;
        idef.il[gmxListedID].iatoms.push_back(parameterIndex);
        idef.il[gmxListedID].iatoms.push_back(index[0]);
        idef.il[gmxListedID].iatoms.push_back(index[1]);
        idef.il[gmxListedID].iatoms.push_back(index[2]);
        idef.il[gmxListedID].iatoms.push_back(index[3]);
    }
}

template<class FiveCenterType>
inline std::enable_if_t<Contains<FiveCenterType, SupportedFiveCenterTypes>{}>
transferIndicesImpl(const ListedTypeData<FiveCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int  parameterIndex = index[5] + offset;
        auto gmxListedID    = FindIndex<FiveCenterType, GmxToNblibMapping>::value;
        idef.il[gmxListedID].iatoms.push_back(parameterIndex);
        idef.il[gmxListedID].iatoms.push_back(index[0]);
        idef.il[gmxListedID].iatoms.push_back(index[1]);
        idef.il[gmxListedID].iatoms.push_back(index[2]);
        idef.il[gmxListedID].iatoms.push_back(index[3]);
        idef.il[gmxListedID].iatoms.push_back(index[4]);
    }
}

template<template<class> class Container, class InteractionType>
inline void transferIndices(const Container<InteractionType>& interactionData,
                            InteractionDefinitions&           idef,
                            [[maybe_unused]] int              offset)
{
    if constexpr (ListedTypeIsImplemented<InteractionType>{})
    {
        transferIndicesImpl(interactionData, idef, offset);
    }
}

} // namespace detail

std::tuple<std::unique_ptr<InteractionDefinitions>, std::unique_ptr<gmx_ffparams_t>>
convertToGmxInteractions(const ListedInteractionData& interactions)
{
    std::unique_ptr<gmx_ffparams_t>         ffparamsHolder = std::make_unique<gmx_ffparams_t>();
    std::unique_ptr<InteractionDefinitions> idefHolder =
            std::make_unique<InteractionDefinitions>(*ffparamsHolder);

    gmx_ffparams_t&         ffparams = *ffparamsHolder;
    InteractionDefinitions& idef     = *idefHolder;

    auto copyParamsOneType = [&ffparams](const auto& interactionElement) {
        detail::transferParameters(interactionElement, ffparams);
    };
    for_each_tuple(copyParamsOneType, interactions);

    // since gmx_ffparams_t.iparams is a flattened vector over all interaction types,
    // we need to compute offsets for each type to know where the parameters for each type start
    // in the flattened iparams vectors
    int                                                       acc = 0;
    std::array<int, std::tuple_size_v<ListedInteractionData>> indexOffsets{ 0 };
    auto extractNIndices = [&indexOffsets, &acc](const auto& interactionElement) {
        constexpr int elementIndex =
                FindIndex<std::decay_t<decltype(interactionElement)>, ListedInteractionData>::value;
        indexOffsets[elementIndex] = acc;
        acc += interactionElement.parameters.size();
    };
    for_each_tuple(extractNIndices, interactions);

    auto copyIndicesOneType = [&idef, &indexOffsets](const auto& interactionElement) {
        constexpr int elementIndex =
                FindIndex<std::decay_t<decltype(interactionElement)>, ListedInteractionData>::value;
        detail::transferIndices(interactionElement, idef, indexOffsets[elementIndex]);
    };
    for_each_tuple(copyIndicesOneType, interactions);

    return std::make_tuple(std::move(idefHolder), std::move(ffparamsHolder));
}


} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CONVERSION_HPP
