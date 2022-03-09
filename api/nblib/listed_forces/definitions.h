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
/*! \inpublicapi \file
 * \brief
 * Definitions for supported nblib listed interaction data, such as bonds, angles, dihedrals, etc
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef NBLIB_LISTEDFORCES_DEFINITIONS_H
#define NBLIB_LISTEDFORCES_DEFINITIONS_H

#include <variant>

#include "nblib/util/traits.hpp"

#include "bondtypes.h"

namespace nblib
{

//***********************************************************************************

/*! \brief These type lists define what interaction types are supported in
 *  -Molecule
 *  -Topology
 *  -ListedForceCalculator
 *
 *  To enable force calculation for your new interaction type that you've added to bondtypes.h,
 *  list your new type here under the appropriate category and make sure that you've added
 *  a kernel in kernels.hpp
 */

using SupportedTwoCenterTypes =
        TypeList<HarmonicBondType, G96BondType, CubicBondType, MorseBondType, FENEBondType, HalfAttractiveQuarticBondType, PairLJType>;
using SupportedThreeCenterTypes =
        TypeList<HarmonicAngle, G96Angle, QuarticAngle, RestrictedAngle, CrossBondBond, CrossBondAngle, LinearAngle>;
using SupportedFourCenterTypes = TypeList<ProperDihedral, ImproperDihedral, RyckaertBellemanDihedral>;
using SupportedFiveCenterTypes = TypeList<Default5Center>;

//***********************************************************************************

using SupportedListedTypes =
        Fuse<SupportedTwoCenterTypes, SupportedThreeCenterTypes, SupportedFourCenterTypes, SupportedFiveCenterTypes>;

//! \brief meta function to map from an Interaction type to the number of interaction centers
template<class Interaction, class = void>
struct NCenter
{
};

//! \brief meta function return value for two-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedTwoCenterTypes>{}>> :
    std::integral_constant<std::size_t, 2>
{
};

//! \brief meta function return value for three-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedThreeCenterTypes>{}>> :
    std::integral_constant<std::size_t, 3>
{
};

//! \brief meta function return value for four-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedFourCenterTypes>{}>> :
    std::integral_constant<std::size_t, 4>
{
};

//! \brief meta function return value for five-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedFiveCenterTypes>{}>> :
    std::integral_constant<std::size_t, 5>
{
};

template<size_t N>
using IndexArray = std::array<int, N>;

/*! \brief encodes the number of integers needed to represent N-center interactions
 *
 *  number of indices to store is the the number of interaction center
 *  plus 1 index for the interaction parameter lookup
 */
template<class Interaction>
using InteractionIndex = IndexArray<NCenter<Interaction>{} + 1>;

// same as InteractionIndex, but just the coordinate indices
template<class Interaction>
using CoordinateIndex = IndexArray<NCenter<Interaction>{}>;

//! \brief container data type for listed interactions
template<class InteractionType>
struct ListedTypeData
{
    using type = InteractionType;

    // vector of unique TwoCenterType instances
    std::vector<InteractionType> parameters;
    // tuple format: <particleID i, particleID j, ..., InteractionInstanceIndex>
    std::vector<InteractionIndex<InteractionType>> indices;
};

using TwoCenterInteraction   = Reduce<std::variant, SupportedTwoCenterTypes>;
using ThreeCenterInteraction = Reduce<std::variant, SupportedThreeCenterTypes>;
using FourCenterInteraction  = Reduce<std::variant, SupportedFourCenterTypes>;
using FiveCenterInteraction  = Reduce<std::variant, SupportedFiveCenterTypes>;

//! This is the complete type that holds all listed interaction data
// result: std::tuple<ListedTypeData<SupportedListedType1>, ...>
using ListedInteractionData = Reduce<std::tuple, Map<ListedTypeData, SupportedListedTypes>>;

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_DEFINITIONS_H
