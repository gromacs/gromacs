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
 * representations of the GMX backed and that of
 * the listed forces parameters for NBLIB
 *
 * To be used to read system topology from TPR files
 * and construct NB-LIB objects from it
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTEDFORCES_CONVERTGMXTONBLIB_H
#define NBLIB_LISTEDFORCES_CONVERTGMXTONBLIB_H

#include <cmath>
#include <cstddef>

#include <algorithm>
#include <array>
#include <memory>
#include <type_traits>
#include <vector>

#include "listed_forces/conversionscommon.h"

#include "gromacs/topology/idef.h"

#include "nblib/listed_forces/bondtypes.h"
#include "nblib/listed_forces/definitions.h"
#include "nblib/particletype.h"
#include "nblib/util/traits.hpp"
#include "nblib/util/util.hpp"

namespace nblib
{

namespace detail
{

template<class InteractionData>
inline void transferParametersGmxToNblib(const t_iparams& /* iparams */,
                                         InteractionData& /* interactionData */)
{
}

//! \brief Harmonic Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                  iparams,
                                         ListedTypeData<HarmonicBondType>& interactionData)
{
    HarmonicBondType harmonicBondType(iparams.harmonic.krA, iparams.harmonic.rA);
    interactionData.parameters.push_back(harmonicBondType);
}

//! \brief G96 Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<G96BondType>& interactionData)
{
    G96BondType g96BondType(iparams.harmonic.krA, std::sqrt(iparams.harmonic.rA));
    interactionData.parameters.push_back(g96BondType);
}

//! \brief FENE Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<FENEBondType>& interactionData)
{
    FENEBondType feneBondType(iparams.fene.kb, iparams.fene.bm);
    interactionData.parameters.push_back(feneBondType);
}

//! \brief Cubic Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<CubicBondType>& interactionData)
{
    CubicBondType cubicBondType(iparams.cubic.kb, iparams.cubic.kcub, iparams.cubic.b0);
    interactionData.parameters.push_back(cubicBondType);
}

//! \brief Morse Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<MorseBondType>& interactionData)
{
    MorseBondType morseBondType(iparams.morse.cbA, iparams.morse.betaA, iparams.morse.b0A);
    interactionData.parameters.push_back(morseBondType);
}

//! \brief Lennard-Jones 1-4 pair interaction parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<PairLJType>& interactionData)
{
    PairLJType pairLjType(C6(iparams.lj14.c6A), C12(iparams.lj14.c12A));
    interactionData.parameters.push_back(pairLjType);
}

//! \brief Harmonic Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<HarmonicAngle>& interactionData)
{
    HarmonicAngle harmonicAngle(iparams.harmonic.krA, Degrees(iparams.harmonic.rA));
    interactionData.parameters.push_back(harmonicAngle);
}

//! \brief G96 Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<G96Angle>& interactionData)
{
    G96Angle g96Angle(iparams.harmonic.krA, Radians(std::acos(iparams.harmonic.rA)));
    interactionData.parameters.push_back(g96Angle);
}

//! \brief Linear Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<LinearAngle>& interactionData)
{
    LinearAngle linearAngle(iparams.linangle.klinA, iparams.linangle.aA);
    interactionData.parameters.push_back(linearAngle);
}

//! \brief Restricted Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                 iparams,
                                         ListedTypeData<RestrictedAngle>& interactionData)
{
    RestrictedAngle restrictedAngle(iparams.harmonic.krA, Degrees(iparams.harmonic.rA));
    interactionData.parameters.push_back(restrictedAngle);
}

//! \brief Quartic Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<QuarticAngle>& interactionData)
{
    QuarticAngle quarticAngle(iparams.qangle.c[0],
                              iparams.qangle.c[1],
                              iparams.qangle.c[2],
                              iparams.qangle.c[3],
                              iparams.qangle.c[4],
                              Degrees(iparams.qangle.theta));
    interactionData.parameters.push_back(quarticAngle);
}

//! \brief Cross Bond-Bond parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<CrossBondBond>& interactionData)
{
    CrossBondBond crossBondBond(iparams.cross_bb.krr, iparams.cross_bb.r1e, iparams.cross_bb.r2e);
    interactionData.parameters.push_back(crossBondBond);
}

//! \brief Cross Bond-Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                iparams,
                                         ListedTypeData<CrossBondAngle>& interactionData)
{
    CrossBondAngle crossBondAngle(
            iparams.cross_ba.krt, iparams.cross_ba.r1e, iparams.cross_ba.r2e, iparams.cross_ba.r3e);
    interactionData.parameters.push_back(crossBondAngle);
}

//! \brief Proper Dihedral parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                iparams,
                                         ListedTypeData<ProperDihedral>& interactionData)
{
    ProperDihedral properDihedral(Degrees(iparams.pdihs.phiA), iparams.pdihs.cpA, iparams.pdihs.mult);
    interactionData.parameters.push_back(properDihedral);
}

//! \brief Ryckaert-Belleman Dihedral parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                          iparams,
                                         ListedTypeData<RyckaertBellemanDihedral>& interactionData)
{
    RyckaertBellemanDihedral ryckaertBellemanDihedral(iparams.rbdihs.rbcA[0],
                                                      iparams.rbdihs.rbcA[1],
                                                      iparams.rbdihs.rbcA[2],
                                                      iparams.rbdihs.rbcA[3],
                                                      iparams.rbdihs.rbcA[4],
                                                      iparams.rbdihs.rbcA[5]);
    interactionData.parameters.push_back(ryckaertBellemanDihedral);
}

template<class TwoCenterType>
inline std::enable_if_t<Contains<TwoCenterType, GmxToNblibMapping>{} && Contains<TwoCenterType, SupportedTwoCenterTypes>{}>
transferIndicesGmxToNblib(const InteractionDefinitions&  idef,
                          ListedTypeData<TwoCenterType>& interactionData,
                          const std::vector<int>&        uniqueParamIndices)
{
    // copy indices
    auto   gmxListedID           = FindIndex<TwoCenterType, GmxToNblibMapping>::value;
    size_t numInteractionsofType = idef.il[gmxListedID].iatoms.size() / 3;
    for (size_t i = 0; i < numInteractionsofType; i++)
    {
        InteractionIndex<TwoCenterType> indices;
        indices[0]    = idef.il[gmxListedID].iatoms[3 * i + 1];
        indices[1]    = idef.il[gmxListedID].iatoms[3 * i + 2];
        auto gmxIndex = idef.il[gmxListedID].iatoms[3 * i];
        auto pointer = std::lower_bound(uniqueParamIndices.begin(), uniqueParamIndices.end(), gmxIndex);
        indices[2] = pointer - uniqueParamIndices.begin();
        interactionData.indices.push_back(indices);
    }
}

template<class ThreeCenterType>
inline std::enable_if_t<Contains<ThreeCenterType, GmxToNblibMapping>{} && Contains<ThreeCenterType, SupportedThreeCenterTypes>{}>
transferIndicesGmxToNblib(const InteractionDefinitions&    idef,
                          ListedTypeData<ThreeCenterType>& interactionData,
                          const std::vector<int>&          uniqueParamIndices)
{
    // copy indices
    auto   gmxListedID           = FindIndex<ThreeCenterType, GmxToNblibMapping>::value;
    size_t numInteractionsofType = idef.il[gmxListedID].iatoms.size() / 4;
    for (size_t i = 0; i < numInteractionsofType; i++)
    {
        InteractionIndex<ThreeCenterType> indices;
        indices[0]    = idef.il[gmxListedID].iatoms[4 * i + 1];
        indices[1]    = idef.il[gmxListedID].iatoms[4 * i + 2];
        indices[2]    = idef.il[gmxListedID].iatoms[4 * i + 3];
        auto gmxIndex = idef.il[gmxListedID].iatoms[4 * i];
        auto pointer = std::lower_bound(uniqueParamIndices.begin(), uniqueParamIndices.end(), gmxIndex);
        indices[3] = pointer - uniqueParamIndices.begin();
        interactionData.indices.push_back(indices);
    }
}

template<class FourCenterType>
inline std::enable_if_t<Contains<FourCenterType, GmxToNblibMapping>{} && Contains<FourCenterType, SupportedFourCenterTypes>{}>
transferIndicesGmxToNblib(const InteractionDefinitions&   idef,
                          ListedTypeData<FourCenterType>& interactionData,
                          const std::vector<int>&         uniqueParamIndices)
{
    // copy indices
    auto   gmxListedID           = FindIndex<FourCenterType, GmxToNblibMapping>::value;
    size_t numInteractionsofType = idef.il[gmxListedID].iatoms.size() / 5;
    for (size_t i = 0; i < numInteractionsofType; i++)
    {
        InteractionIndex<FourCenterType> indices;
        indices[0]    = idef.il[gmxListedID].iatoms[5 * i + 1];
        indices[1]    = idef.il[gmxListedID].iatoms[5 * i + 2];
        indices[2]    = idef.il[gmxListedID].iatoms[5 * i + 3];
        indices[3]    = idef.il[gmxListedID].iatoms[5 * i + 4];
        auto gmxIndex = idef.il[gmxListedID].iatoms[5 * i];
        auto pointer = std::lower_bound(uniqueParamIndices.begin(), uniqueParamIndices.end(), gmxIndex);
        indices[4] = pointer - uniqueParamIndices.begin();
        interactionData.indices.push_back(indices);
    }
}

// TODO: 5 center interactions like CMAP params are stored differently
//       Concrete logic to be implemented when supported in NBLIB
template<class FiveCenterType>
inline std::enable_if_t<Contains<FiveCenterType, GmxToNblibMapping>{} && Contains<FiveCenterType, SupportedFiveCenterTypes>{}>
transferIndicesGmxToNblib(const InteractionDefinitions& /* idef */,
                          ListedTypeData<FiveCenterType>& /* interactionData */,
                          const std::vector<int>& /* uniqueParamIndices */)
{
}

} // namespace detail

ListedInteractionData convertToNblibInteractions(const InteractionDefinitions& interactionDefinitions)
{
    ListedInteractionData interactions;

    auto transferParamsAndIndices = [&interactionDefinitions](auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;

        if constexpr (Contains<InteractionType, GmxToNblibMapping>{})
        {
            constexpr size_t typeID     = FindIndex<InteractionType, GmxToNblibMapping>::value;
            int              numCenters = NCenter<InteractionType>{};
            auto numInstancesOfType = interactionDefinitions.il[typeID].size() / (numCenters + 1);
            // calculate number of unique parameter sets
            // initiate with a vector of param set IDs
            std::vector<int> uniqueParamIndices(numInstancesOfType);
            for (auto i = 0; i < numInstancesOfType; i++)
            {
                uniqueParamIndices[i] = interactionDefinitions.il[typeID].iatoms[i * (numCenters + 1)];
            }
            // remove duplicate param IDs
            std::sort(uniqueParamIndices.begin(), uniqueParamIndices.end());
            uniqueParamIndices.erase(std::unique(uniqueParamIndices.begin(), uniqueParamIndices.end()),
                                     uniqueParamIndices.end());

            size_t numUniqueParamSets = uniqueParamIndices.size();
            // start filling params and indices into interactions
            // copy param sets
            for (size_t i = 0; i < numUniqueParamSets; i++)
            {
                detail::transferParametersGmxToNblib(
                        interactionDefinitions.iparams[uniqueParamIndices[i]], interactionElement);
            }
            // copy index sets
            detail::transferIndicesGmxToNblib(interactionDefinitions, interactionElement, uniqueParamIndices);
        }
    };
    for_each_tuple(transferParamsAndIndices, interactions);

    return interactions;
}

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CONVERTGMXTONBLIB_H
