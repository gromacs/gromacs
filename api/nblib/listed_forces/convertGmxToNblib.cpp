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
 * and construct NBLIB objects from it
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTEDFORCES_CONVERTGMXTONBLIB_H
#define NBLIB_LISTEDFORCES_CONVERTGMXTONBLIB_H

#include <memory>

#include "nblib/listed_forces/conversionscommon.h"

namespace nblib
{

//! \brief meta function return value for Urey-Bradley interactions
template<>
struct NCenter<UreyBradley> : std::integral_constant<std::size_t, 3>
{
};

namespace detail
{

template<class InteractionData>
inline void transferParametersGmxToNblib(const t_iparams& /* iparams */,
                                         InteractionData& /* interactionData */)
{
}

//! \brief Position restraints parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                    iparams,
                                         ListedTypeData<PositionRestraints>& interactionData)
{
    PositionRestraints positionRestraintsA(iparams.posres.pos0A[0],
                                           iparams.posres.pos0A[1],
                                           iparams.posres.pos0A[2],
                                           iparams.posres.fcA[0],
                                           iparams.posres.fcA[1],
                                           iparams.posres.fcA[2]);
    interactionData.parametersA.push_back(positionRestraintsA);

    PositionRestraints positionRestraintsB(iparams.posres.pos0B[0],
                                           iparams.posres.pos0B[1],
                                           iparams.posres.pos0B[2],
                                           iparams.posres.fcB[0],
                                           iparams.posres.fcB[1],
                                           iparams.posres.fcB[2]);
    interactionData.parametersB.push_back(positionRestraintsB);
}

//! \brief Harmonic Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                  iparams,
                                         ListedTypeData<HarmonicBondType>& interactionData)
{
    HarmonicBondType harmonicBondTypeA(iparams.harmonic.krA, iparams.harmonic.rA);
    interactionData.parametersA.push_back(harmonicBondTypeA);

    HarmonicBondType harmonicBondTypeB(iparams.harmonic.krB, iparams.harmonic.rB);
    interactionData.parametersB.push_back(harmonicBondTypeB);
}

//! \brief G96 Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<G96BondType>& interactionData)
{
    G96BondType g96BondTypeA(iparams.harmonic.krA, std::sqrt(iparams.harmonic.rA));
    interactionData.parametersA.push_back(g96BondTypeA);

    G96BondType g96BondTypeB(iparams.harmonic.krB, std::sqrt(iparams.harmonic.rB));
    interactionData.parametersB.push_back(g96BondTypeB);
}

//! \brief FENE Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<FENEBondType>& interactionData)
{
    FENEBondType feneBondType(iparams.fene.kb, iparams.fene.bm);
    interactionData.parametersA.push_back(feneBondType);
    interactionData.parametersB.push_back(feneBondType);
}

//! \brief Cubic Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<CubicBondType>& interactionData)
{
    CubicBondType cubicBondType(iparams.cubic.kb, iparams.cubic.kcub, iparams.cubic.b0);
    interactionData.parametersA.push_back(cubicBondType);
    interactionData.parametersB.push_back(cubicBondType);
}

//! \brief Morse Bonds parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<MorseBondType>& interactionData)
{
    MorseBondType morseBondTypeA(iparams.morse.cbA, iparams.morse.betaA, iparams.morse.b0A);
    interactionData.parametersA.push_back(morseBondTypeA);

    MorseBondType morseBondTypeB(iparams.morse.cbB, iparams.morse.betaB, iparams.morse.b0B);
    interactionData.parametersB.push_back(morseBondTypeB);
}

//! \brief Lennard-Jones 1-4 pair interaction parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<PairLJType>& interactionData)
{
    PairLJType pairLjTypeA(C6(iparams.lj14.c6A), C12(iparams.lj14.c12A));
    interactionData.parametersA.push_back(pairLjTypeA);

    PairLJType pairLjTypeB(C6(iparams.lj14.c6B), C12(iparams.lj14.c12B));
    interactionData.parametersB.push_back(pairLjTypeB);
}

//! \brief Lennard-Jones 1-4 pair interaction parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                  iparams,
                                         ListedTypeData<PairLJChargeType>& interactionData)
{
    PairLJChargeType pairLjChargeType(C6(iparams.ljc14.c6),
                                      C12(iparams.ljc14.c12),
                                      iparams.ljc14.qi,
                                      iparams.ljc14.qj,
                                      iparams.ljc14.fqq);
    interactionData.parametersA.push_back(pairLjChargeType);
    interactionData.parametersB.push_back(pairLjChargeType);
}

//! \brief Simple polarization parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                    iparams,
                                         ListedTypeData<SimplePolarization>& interactionData)
{
    SimplePolarization simplePolarization(iparams.polarize.alpha);
    interactionData.parametersA.push_back(simplePolarization);
    interactionData.parametersB.push_back(simplePolarization);
}

//! \brief Harmonic Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<HarmonicAngle>& interactionData)
{
    HarmonicAngle harmonicAngleA(iparams.harmonic.krA, Degrees(iparams.harmonic.rA));
    interactionData.parametersA.push_back(harmonicAngleA);

    HarmonicAngle harmonicAngleB(iparams.harmonic.krB, Degrees(iparams.harmonic.rB));
    interactionData.parametersB.push_back(harmonicAngleB);
}

//! \brief G96 Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<G96Angle>& interactionData)
{
    G96Angle g96AngleA(iparams.harmonic.krA, Radians(std::acos(iparams.harmonic.rA)));
    interactionData.parametersA.push_back(g96AngleA);

    G96Angle g96AngleB(iparams.harmonic.krB, Radians(std::acos(iparams.harmonic.rB)));
    interactionData.parametersB.push_back(g96AngleB);
}

//! \brief Linear Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<LinearAngle>& interactionData)
{
    LinearAngle linearAngleA(iparams.linangle.klinA, iparams.linangle.aA);
    interactionData.parametersA.push_back(linearAngleA);

    LinearAngle linearAngleB(iparams.linangle.klinB, iparams.linangle.aB);
    interactionData.parametersB.push_back(linearAngleB);
}

//! \brief Restricted Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                 iparams,
                                         ListedTypeData<RestrictedAngle>& interactionData)
{
    RestrictedAngle restrictedAngleA(iparams.harmonic.krA, Degrees(iparams.harmonic.rA));
    interactionData.parametersA.push_back(restrictedAngleA);

    RestrictedAngle restrictedAngleB(iparams.harmonic.krB, Degrees(iparams.harmonic.rB));
    interactionData.parametersB.push_back(restrictedAngleB);
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
    interactionData.parametersA.push_back(quarticAngle);
    interactionData.parametersB.push_back(quarticAngle);
}

//! \brief Cross Bond-Bond parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams& iparams, ListedTypeData<CrossBondBond>& interactionData)
{
    CrossBondBond crossBondBond(iparams.cross_bb.krr, iparams.cross_bb.r1e, iparams.cross_bb.r2e);
    interactionData.parametersA.push_back(crossBondBond);
    interactionData.parametersB.push_back(crossBondBond);
}

//! \brief Cross Bond-Angle parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                iparams,
                                         ListedTypeData<CrossBondAngle>& interactionData)
{
    CrossBondAngle crossBondAngle(
            iparams.cross_ba.krt, iparams.cross_ba.r1e, iparams.cross_ba.r2e, iparams.cross_ba.r3e);
    interactionData.parametersA.push_back(crossBondAngle);
    interactionData.parametersB.push_back(crossBondAngle);
}

//! \brief Proper Dihedral parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                iparams,
                                         ListedTypeData<ProperDihedral>& interactionData)
{
    ProperDihedral properDihedralA(Degrees(iparams.pdihs.phiA), iparams.pdihs.cpA, iparams.pdihs.mult);
    interactionData.parametersA.push_back(properDihedralA);

    ProperDihedral properDihedralB(Degrees(iparams.pdihs.phiB), iparams.pdihs.cpB, iparams.pdihs.mult);
    interactionData.parametersB.push_back(properDihedralA);
}

//! \brief Improper proper Dihedral parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                        iparams,
                                         ListedTypeData<ImproperProperDihedral>& interactionData)
{
    ImproperProperDihedral properDihedralA(
            Degrees(iparams.pdihs.phiA), iparams.pdihs.cpA, iparams.pdihs.mult);
    interactionData.parametersA.push_back(properDihedralA);

    ImproperProperDihedral properDihedralB(
            Degrees(iparams.pdihs.phiB), iparams.pdihs.cpB, iparams.pdihs.mult);
    interactionData.parametersB.push_back(properDihedralA);
}

//! \brief Improper Dihedral parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                  iparams,
                                         ListedTypeData<ImproperDihedral>& interactionData)
{
    ImproperDihedral improperDihedralA(Degrees(iparams.harmonic.rA), iparams.harmonic.krA);
    interactionData.parametersA.push_back(improperDihedralA);

    ImproperDihedral improperDihedralB(Degrees(iparams.harmonic.rB), iparams.harmonic.krB);
    interactionData.parametersB.push_back(improperDihedralA);
}

//! \brief Ryckaert-Bellman Dihedral parameter transfer function
template<>
inline void transferParametersGmxToNblib(const t_iparams&                          iparams,
                                         ListedTypeData<RyckaertBellemanDihedral>& interactionData)
{
    RyckaertBellemanDihedral ryckaertBellemanDihedralA(iparams.rbdihs.rbcA[0],
                                                       iparams.rbdihs.rbcA[1],
                                                       iparams.rbdihs.rbcA[2],
                                                       iparams.rbdihs.rbcA[3],
                                                       iparams.rbdihs.rbcA[4],
                                                       iparams.rbdihs.rbcA[5]);
    interactionData.parametersA.push_back(ryckaertBellemanDihedralA);

    RyckaertBellemanDihedral ryckaertBellemanDihedralB(iparams.rbdihs.rbcB[0],
                                                       iparams.rbdihs.rbcB[1],
                                                       iparams.rbdihs.rbcB[2],
                                                       iparams.rbdihs.rbcB[3],
                                                       iparams.rbdihs.rbcB[4],
                                                       iparams.rbdihs.rbcB[5]);
    interactionData.parametersB.push_back(ryckaertBellemanDihedralB);
}

template<class OneCenterType>
inline std::enable_if_t<Contains<OneCenterType, GmxToNblibMapping>{} && NCenter<OneCenterType>{} == 1>
transferIndicesGmxToNblib(const InteractionDefinitions&  idef,
                          ListedTypeData<OneCenterType>& interactionData,
                          const std::vector<int>&        uniqueParamIndices)
{
    // copy indices
    auto   gmxListedID           = FindIndex<OneCenterType, GmxToNblibMapping>::value;
    size_t numInteractionsofType = idef.il[gmxListedID].iatoms.size() / 2;
    for (size_t i = 0; i < numInteractionsofType; i++)
    {
        InteractionIndex<OneCenterType> indices;
        indices[0]    = idef.il[gmxListedID].iatoms[2 * i + 1];
        auto gmxIndex = idef.il[gmxListedID].iatoms[2 * i];
        auto pointer = std::lower_bound(uniqueParamIndices.begin(), uniqueParamIndices.end(), gmxIndex);
        indices[1] = pointer - uniqueParamIndices.begin();
        interactionData.indices.push_back(indices);
    }
}

template<class TwoCenterType>
inline std::enable_if_t<Contains<TwoCenterType, GmxToNblibMapping>{} && NCenter<TwoCenterType>{} == 2>
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
inline std::enable_if_t<Contains<ThreeCenterType, GmxToNblibMapping>{} && NCenter<ThreeCenterType>{} == 3>
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
inline std::enable_if_t<Contains<FourCenterType, GmxToNblibMapping>{} && NCenter<FourCenterType>{} == 4>
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
inline std::enable_if_t<Contains<FiveCenterType, GmxToNblibMapping>{} && NCenter<FiveCenterType>{} == 5>
transferIndicesGmxToNblib(const InteractionDefinitions& /* idef */,
                          ListedTypeData<FiveCenterType>& /* interactionData */,
                          const std::vector<int>& /* uniqueParamIndices */)
{
}

template<class InteractionType>
inline std::vector<int> calculateUniqueParamIndices(const InteractionDefinitions& interactionDefinitions)
{
    constexpr size_t typeID             = FindIndex<InteractionType, GmxToNblibMapping>::value;
    auto             indexVector        = interactionDefinitions.il[typeID].iatoms;
    int              numCenters         = NCenter<InteractionType>{};
    auto             numInstancesOfType = indexVector.size() / (numCenters + 1);

    // initiate with a vector of param set IDs
    std::vector<int> uniqueParamIndices(numInstancesOfType);
    for (size_t i = 0; i < numInstancesOfType; i++)
    {
        uniqueParamIndices[i] = indexVector[i * (numCenters + 1)];
    }

    // remove duplicate param IDs
    std::sort(uniqueParamIndices.begin(), uniqueParamIndices.end());
    uniqueParamIndices.erase(std::unique(uniqueParamIndices.begin(), uniqueParamIndices.end()),
                             uniqueParamIndices.end());

    return uniqueParamIndices;
}

//! \brief Transfer UB interactions from GMX as a separate Harmonic Bond and Harmonic Angle in NBLIB
static void transferUreyBradley(const InteractionDefinitions& interactionDefinitions,
                                ListedInteractionData&        interactions)
{
    // check if UB interactions exist in the system
    constexpr size_t gmxListedID = FindIndex<UreyBradley, GmxToNblibMapping>::value;
    if (!interactionDefinitions.il[gmxListedID].iatoms.empty())
    {
        auto uniqueParamIndices = detail::calculateUniqueParamIndices<UreyBradley>(interactionDefinitions);
        auto numHarmonicBonds  = pickType<HarmonicBondType>(interactions).parametersA.size();
        auto numHarmonicAngles = pickType<HarmonicAngle>(interactions).parametersA.size();

        // copy param sets
        for (size_t i = 0; i < uniqueParamIndices.size(); i++)
        {
            auto param = interactionDefinitions.iparams[uniqueParamIndices[i]];
            // push harmonic bond
            HarmonicBondType harmonicBondTypeA(param.u_b.kUBA, param.u_b.r13A);
            pickType<HarmonicBondType>(interactions).parametersA.push_back(harmonicBondTypeA);
            HarmonicBondType harmonicBondTypeB(param.u_b.kUBB, param.u_b.r13B);
            pickType<HarmonicBondType>(interactions).parametersB.push_back(harmonicBondTypeB);

            // push harmonic angle
            HarmonicAngle harmonicAngleA(param.u_b.kthetaA, Degrees(param.u_b.thetaA));
            pickType<HarmonicAngle>(interactions).parametersA.push_back(harmonicAngleA);
            HarmonicAngle harmonicAngleB(param.u_b.kthetaB, Degrees(param.u_b.thetaB));
            pickType<HarmonicAngle>(interactions).parametersB.push_back(harmonicAngleB);
        }

        // copy index sets
        size_t numUreyBradley = interactionDefinitions.il[gmxListedID].iatoms.size() / 4;
        for (size_t i = 0; i < numUreyBradley; i++)
        {
            InteractionIndex<HarmonicBondType> hbIndices;
            InteractionIndex<HarmonicAngle>    haIndices;

            auto gmxIndex = interactionDefinitions.il[gmxListedID].iatoms[4 * i];
            auto it = std::lower_bound(uniqueParamIndices.begin(), uniqueParamIndices.end(), gmxIndex);

            hbIndices[0] = interactionDefinitions.il[gmxListedID].iatoms[4 * i + 1];
            hbIndices[1] = interactionDefinitions.il[gmxListedID].iatoms[4 * i + 3];
            hbIndices[2] = it - uniqueParamIndices.begin() + numHarmonicBonds;

            haIndices[0] = interactionDefinitions.il[gmxListedID].iatoms[4 * i + 1];
            haIndices[1] = interactionDefinitions.il[gmxListedID].iatoms[4 * i + 2];
            haIndices[2] = interactionDefinitions.il[gmxListedID].iatoms[4 * i + 3];
            haIndices[3] = it - uniqueParamIndices.begin() + numHarmonicAngles;

            pickType<HarmonicBondType>(interactions).indices.push_back(hbIndices);
            pickType<HarmonicAngle>(interactions).indices.push_back(haIndices);
        }
    }
}

} // namespace detail

ListedInteractionData convertToNblibInteractions(const InteractionDefinitions& interactionDefinitions)
{
    ListedInteractionData interactions;

    auto transferParamsAndIndices = [&interactionDefinitions](auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;

        if constexpr (Contains<InteractionType, GmxToNblibMapping>{})
        {
            // calculate vector of unique parameter sets
            auto uniqueParamIndices =
                    detail::calculateUniqueParamIndices<InteractionType>(interactionDefinitions);

            const auto& parameterSource = (Contains<InteractionType, RestraintTypes>{})
                                                  ? interactionDefinitions.iparams_posres
                                                  : interactionDefinitions.iparams;
            // copy param sets
            for (size_t i = 0; i < uniqueParamIndices.size(); i++)
            {
                detail::transferParametersGmxToNblib(parameterSource[uniqueParamIndices[i]],
                                                     interactionElement);
            }
            // copy index sets
            detail::transferIndicesGmxToNblib(interactionDefinitions, interactionElement, uniqueParamIndices);
        }
    };
    for_each_tuple(transferParamsAndIndices, interactions);

    detail::transferUreyBradley(interactionDefinitions, interactions);

    return interactions;
}

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CONVERTGMXTONBLIB_H
