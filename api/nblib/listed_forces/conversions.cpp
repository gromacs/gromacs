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

#include <memory>
#include <vector>

#include "gromacs/topology/forcefieldparameters.h"

#include "nblib/listed_forces/conversionscommon.h"

namespace nblib
{

namespace detail
{

//! \brief fall-back for nblib types that do not exist in GROMACS
template<class InteractionData>
inline void transferParameters(const InteractionData&, std::vector<t_iparams>&)
{
}

//! \brief Position restraints parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<PositionRestraints>& interactions,
                               std::vector<t_iparams>&                   iparams_posres)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Position restraints size mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        for (int j = 0; j < DIM; j++)
        {
            param.posres.fcA[j]   = interactions.parametersA[i].forceConstant(j);
            param.posres.pos0A[j] = interactions.parametersA[i].position0(j);
        }
        if (!interactions.parametersB.empty())
        {
            for (int j = 0; j < DIM; j++)
            {
                param.posres.fcB[j]   = interactions.parametersB[i].forceConstant(j);
                param.posres.pos0B[j] = interactions.parametersB[i].position0(j);
            }
        }
        else
        {
            for (int j = 0; j < DIM; j++)
            {
                param.posres.fcB[j]   = 0;
                param.posres.pos0B[j] = 0;
            }
        }
        iparams_posres.push_back(param);
    }
}

//! \brief Harmonic bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<HarmonicBondType>& interactions,
                               std::vector<t_iparams>&                 gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Harmonic bond interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactions.parametersA[i].forceConstant();
        param.harmonic.rA  = interactions.parametersA[i].equilConstant();
        if (!interactions.parametersB.empty())
        {
            param.harmonic.krB = interactions.parametersB[i].forceConstant();
            param.harmonic.rB  = interactions.parametersB[i].equilConstant();
        }
        else
        {
            param.harmonic.krB = 0.0;
            param.harmonic.rB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief G96 bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<G96BondType>& interactions,
                               std::vector<t_iparams>&            gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("G96 bond interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactions.parametersA[i].forceConstant();
        param.harmonic.rA  = interactions.parametersA[i].equilConstant();
        if (!interactions.parametersB.empty())
        {
            param.harmonic.krB = interactions.parametersB[i].forceConstant();
            param.harmonic.rB  = interactions.parametersB[i].equilConstant();
        }
        else
        {
            param.harmonic.krB = 0.0;
            param.harmonic.rB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief FENE bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<FENEBondType>& interactions,
                               std::vector<t_iparams>&             gmx_params)
{
    for (const auto& bond : interactions.parametersA)
    {
        t_iparams param;
        param.fene.kb = bond.forceConstant();
        param.fene.bm = bond.equilConstant();
        gmx_params.push_back(param);
    }
}

//! \brief Cubic bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<CubicBondType>& interactions,
                               std::vector<t_iparams>&              gmx_params)
{
    for (const auto& bond : interactions.parametersA)
    {
        t_iparams param;
        param.cubic.kcub = bond.cubicForceConstant();
        param.cubic.kb   = bond.quadraticForceConstant();
        param.cubic.b0   = bond.equilDistance();
        gmx_params.push_back(param);
    }
}

//! \brief Morse bond parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<MorseBondType>& interactions,
                               std::vector<t_iparams>&              gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Morse bond interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.morse.b0A   = interactions.parametersA[i].equilDistance();
        param.morse.betaA = interactions.parametersA[i].exponent();
        param.morse.cbA   = interactions.parametersA[i].forceConstant();
        if (!interactions.parametersB.empty())
        {
            param.morse.b0B   = interactions.parametersB[i].equilDistance();
            param.morse.betaB = interactions.parametersB[i].exponent();
            param.morse.cbB   = interactions.parametersB[i].forceConstant();
        }
        else
        {
            param.morse.b0B   = 0.0;
            param.morse.betaB = 0.0;
            param.morse.cbB   = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Pair LJ 1-4 interaction parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<PairLJType>& interactions,
                               std::vector<t_iparams>&           gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("LJ1-4 pair interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.lj14.c6A  = interactions.parametersA[i].c6();
        param.lj14.c12A = interactions.parametersA[i].c12();
        if (!interactions.parametersB.empty())
        {
            param.lj14.c6B  = interactions.parametersB[i].c6();
            param.lj14.c12B = interactions.parametersB[i].c12();
        }
        else
        {
            param.lj14.c6B  = 0.0;
            param.lj14.c12B = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Charged pair LJ 1-4 interaction parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<PairLJChargeType>& interactions,
                               std::vector<t_iparams>&                 gmx_params)
{
    for (const auto& pair : interactions.parametersA)
    {
        t_iparams param;
        param.ljc14.c6  = pair.c6();
        param.ljc14.c12 = pair.c12();
        param.ljc14.qi  = pair.qi();
        param.ljc14.qj  = pair.qj();
        param.ljc14.fqq = pair.ff();
        gmx_params.push_back(param);
    }
}

//! \brief Simple polarization parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<SimplePolarization>& interactions,
                               std::vector<t_iparams>&                   gmx_params)
{
    for (const auto& bond : interactions.parametersA)
    {
        t_iparams param;
        param.polarize.alpha = bond.alpha();
        gmx_params.push_back(param);
    }
}

//! \brief Harmonic angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<HarmonicAngle>& interactions,
                               std::vector<t_iparams>&              gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Harmonic angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactions.parametersA[i].forceConstant();
        param.harmonic.rA  = interactions.parametersA[i].equilConstant() / DEG2RAD;
        if (!interactions.parametersB.empty())
        {
            param.harmonic.krB = interactions.parametersB[i].forceConstant();
            param.harmonic.rB  = interactions.parametersB[i].equilConstant() / DEG2RAD;
        }
        else
        {
            param.harmonic.krB = 0.0;
            param.harmonic.rB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief G96 angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<G96Angle>& interactions, std::vector<t_iparams>& gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("G96 angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactions.parametersA[i].forceConstant();
        param.harmonic.rA  = interactions.parametersA[i].equilConstant();
        if (!interactions.parametersB.empty())
        {
            param.harmonic.krB = interactions.parametersB[i].forceConstant();
            param.harmonic.rB  = interactions.parametersB[i].equilConstant();
        }
        else
        {
            param.harmonic.krB = 0.0;
            param.harmonic.rB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Linear angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<LinearAngle>& interactions,
                               std::vector<t_iparams>&            gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Linear angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.linangle.klinA = interactions.parametersA[i].forceConstant();
        param.linangle.aA    = interactions.parametersA[i].equilConstant();
        if (!interactions.parametersB.empty())
        {
            param.linangle.klinB = interactions.parametersB[i].forceConstant();
            param.linangle.aB    = interactions.parametersB[i].equilConstant();
        }
        else
        {
            param.linangle.klinB = 0.0;
            param.linangle.aB    = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Restricted angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<RestrictedAngle>& interactions,
                               std::vector<t_iparams>&                gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Restricted angle interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactions.parametersA[i].forceConstant();
        param.harmonic.rA  = (std::acos(interactions.parametersA[i].equilConstant())) / DEG2RAD;
        if (!interactions.parametersB.empty())
        {
            param.harmonic.krB = interactions.parametersB[i].forceConstant();
            param.harmonic.rB  = (std::acos(interactions.parametersB[i].equilConstant())) / DEG2RAD;
        }
        else
        {
            param.harmonic.krB = 0.0;
            param.harmonic.rB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Quartic angle parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<QuarticAngle>& interactions,
                               std::vector<t_iparams>&             gmx_params)
{
    for (const auto& angle : interactions.parametersA)
    {
        t_iparams param;
        param.qangle.theta = angle.equilConstant() / DEG2RAD;
        param.qangle.c[0]  = angle.forceConstant(0);
        param.qangle.c[1]  = angle.forceConstant(1);
        param.qangle.c[2]  = angle.forceConstant(2);
        param.qangle.c[3]  = angle.forceConstant(3);
        param.qangle.c[4]  = angle.forceConstant(4);
        gmx_params.push_back(param);
    }
}

//! \brief Cross Bond-Bond interaction parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<CrossBondBond>& interactions,
                               std::vector<t_iparams>&              gmx_params)
{
    for (const auto& bond : interactions.parametersA)
    {
        t_iparams param;
        param.cross_bb.krr = bond.forceConstant();
        param.cross_bb.r1e = bond.equilDistanceIJ();
        param.cross_bb.r2e = bond.equilDistanceKJ();
        gmx_params.push_back(param);
    }
}

//! \brief Cross Bond-Angle interaction parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<CrossBondAngle>& interactions,
                               std::vector<t_iparams>&               gmx_params)
{
    for (const auto& bond : interactions.parametersA)
    {
        t_iparams param;
        param.cross_ba.krt = bond.forceConstant();
        param.cross_ba.r1e = bond.equilDistanceIJ();
        param.cross_ba.r2e = bond.equilDistanceKJ();
        param.cross_ba.r3e = bond.equilDistanceIK();
        gmx_params.push_back(param);
    }
}

//! \brief Proper dihedral parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<ProperDihedral>& interactions,
                               std::vector<t_iparams>&               gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Proper dihedral interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.pdihs.phiA = interactions.parametersA[i].equilDistance() / DEG2RAD;
        param.pdihs.cpA  = interactions.parametersA[i].forceConstant();
        param.pdihs.mult = interactions.parametersA[i].multiplicity();
        if (!interactions.parametersB.empty())
        {
            param.pdihs.phiB = interactions.parametersB[i].equilDistance() / DEG2RAD;
            param.pdihs.cpB  = interactions.parametersB[i].forceConstant();
        }
        else
        {
            param.pdihs.phiB = 0.0;
            param.pdihs.cpB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Proper dihedral parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<ImproperProperDihedral>& interactions,
                               std::vector<t_iparams>&                       gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Proper dihedral interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.pdihs.phiA = interactions.parametersA[i].equilDistance() / DEG2RAD;
        param.pdihs.cpA  = interactions.parametersA[i].forceConstant();
        param.pdihs.mult = interactions.parametersA[i].multiplicity();
        if (!interactions.parametersB.empty())
        {
            param.pdihs.phiB = interactions.parametersB[i].equilDistance() / DEG2RAD;
            param.pdihs.cpB  = interactions.parametersB[i].forceConstant();
        }
        else
        {
            param.pdihs.phiB = 0.0;
            param.pdihs.cpB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Improper dihedral parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<ImproperDihedral>& interactions,
                               std::vector<t_iparams>&                 gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Improper dihedral interactions array mismatch for A & B");
    }
    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.harmonic.krA = interactions.parametersA[i].forceConstant();
        param.harmonic.rA  = interactions.parametersA[i].equilConstant() / DEG2RAD;
        if (!interactions.parametersB.empty())
        {
            param.harmonic.krB = interactions.parametersB[i].forceConstant();
            param.harmonic.rB  = interactions.parametersB[i].equilConstant() / DEG2RAD;
        }
        else
        {
            param.harmonic.krB = 0.0;
            param.harmonic.rB  = 0.0;
        }
        gmx_params.push_back(param);
    }
}

//! \brief Ryckaert-Belleman dihedral parameter conversion function
template<>
inline void transferParameters(const ListedTypeData<RyckaertBellemanDihedral>& interactions,
                               std::vector<t_iparams>&                         gmx_params)
{
    if (interactions.parametersA.size() != interactions.parametersB.size()
        && !interactions.parametersB.empty())
    {
        throw InputException("Ryckaert-Belleman dihedral interactions array mismatch for A & B");
    }

    for (size_t i = 0; i < interactions.parametersA.size(); i++)
    {
        t_iparams param;
        param.rbdihs.rbcA[0] = interactions.parametersA[i][0];
        param.rbdihs.rbcA[1] = interactions.parametersA[i][1];
        param.rbdihs.rbcA[2] = interactions.parametersA[i][2];
        param.rbdihs.rbcA[3] = interactions.parametersA[i][3];
        param.rbdihs.rbcA[4] = interactions.parametersA[i][4];
        param.rbdihs.rbcA[5] = interactions.parametersA[i][5];
        if (!interactions.parametersB.empty())
        {
            param.rbdihs.rbcB[0] = interactions.parametersB[i][0];
            param.rbdihs.rbcB[1] = interactions.parametersB[i][1];
            param.rbdihs.rbcB[2] = interactions.parametersB[i][2];
            param.rbdihs.rbcB[3] = interactions.parametersB[i][3];
            param.rbdihs.rbcB[4] = interactions.parametersB[i][4];
            param.rbdihs.rbcB[5] = interactions.parametersB[i][5];
        }
        else
        {
            param.rbdihs.rbcB[0] = 0.0;
            param.rbdihs.rbcB[1] = 0.0;
            param.rbdihs.rbcB[2] = 0.0;
            param.rbdihs.rbcB[3] = 0.0;
            param.rbdihs.rbcB[4] = 0.0;
            param.rbdihs.rbcB[5] = 0.0;
        }
        gmx_params.push_back(param);
    }
}

template<class OneCenterType>
inline std::enable_if_t<NCenter<OneCenterType>{} == 1>
transferIndicesImpl(const ListedTypeData<OneCenterType>& interactions, InteractionDefinitions& idef, int offset)
{
    for (const auto& index : interactions.indices)
    {
        int  parameterIndex = index[1] + offset;
        auto gmxListedID    = FindIndex<OneCenterType, GmxToNblibMapping>::value;
        idef.il[gmxListedID].iatoms.push_back(parameterIndex);
        idef.il[gmxListedID].iatoms.push_back(index[0]);
    }
}

template<class TwoCenterType>
inline std::enable_if_t<NCenter<TwoCenterType>{} == 2>
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
inline std::enable_if_t<NCenter<ThreeCenterType>{} == 3>
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
inline std::enable_if_t<NCenter<FourCenterType>{} == 4>
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
inline std::enable_if_t<NCenter<FiveCenterType>{} == 5>
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

    auto copyParamsOneType = [&ffparams, &idef](const auto& interactionElement) {
        using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;

        constexpr bool isRestraint = Contains<InteractionType, RestraintTypes>{};
        // Necessary to specially copy the position restrains data to a specific idef member
        auto& gmxParams = isRestraint ? idef.iparams_posres : ffparams.iparams;
        detail::transferParameters(interactionElement, gmxParams);

        int numParameters   = isRestraint ? 0 : interactionElement.parametersA.size();
        int newFuncTypeSize = ffparams.functype.size() + numParameters;
        ffparams.functype.resize(newFuncTypeSize);
        constexpr int gmxTypeId = FindIndex<InteractionType, GmxToNblibMapping>{};
        std::fill(ffparams.functype.end() - numParameters, ffparams.functype.end(), gmxTypeId);
    };
    for_each_tuple(copyParamsOneType, interactions);

    // since gmx_ffparams_t.iparams is a flattened vector over all interaction types,
    // we need to compute offsets for each type to know where the parameters for each type start
    // in the flattened iparams vectors
    int paramOffset = 0;
    // position restraints have their own parameter container in idef, therefore we need a separate offset
    int posresOffset = 0;

    std::array<int, std::tuple_size_v<ListedInteractionData>> indexOffsets{ 0 };
    auto extractNIndicesBase = [&indexOffsets, &paramOffset, &posresOffset](const auto& interactionElement) {
        constexpr int elementIndex =
                FindIndex<std::decay_t<decltype(interactionElement)>, ListedInteractionData>::value;
        using InteractionType = TypeListElement_t<elementIndex, AllListedTypes>;

        if constexpr (Contains<InteractionType, RestraintTypes>{})
        {
            indexOffsets[elementIndex] = posresOffset;
            posresOffset += interactionElement.parametersA.size();
        }
        else
        {
            indexOffsets[elementIndex] = paramOffset;
            paramOffset += interactionElement.parametersA.size();
        }
    };
    for_each_tuple(extractNIndicesBase, interactions);

    auto copyIndicesOneTypeBase = [&idef, &indexOffsets](const auto& interactionElement) {
        constexpr int elementIndex =
                FindIndex<std::decay_t<decltype(interactionElement)>, ListedInteractionData>::value;
        detail::transferIndices(interactionElement, idef, indexOffsets[elementIndex]);
    };
    for_each_tuple(copyIndicesOneTypeBase, interactions);

    return std::make_tuple(std::move(idefHolder), std::move(ffparamsHolder));
}


} // namespace nblib

#endif // NBLIB_LISTEDFORCES_CONVERSION_HPP
