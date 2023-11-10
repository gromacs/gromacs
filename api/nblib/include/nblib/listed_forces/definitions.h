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

#include <type_traits>
#include <variant>
#include <vector>

#include "nblib/listed_forces/bondtypes.h"
#include "nblib/util/array.hpp"
#include "nblib/util/traits.hpp"

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

//! \brief basic listed types without charges
using SupportedTwoCenterTypes =
        TypeList<HarmonicBondType, G96BondType, CubicBondType, MorseBondType, FENEBondType, HalfAttractiveQuarticBondType>;
using SupportedThreeCenterTypes =
        TypeList<HarmonicAngle, G96Angle, QuarticAngle, RestrictedAngle, CrossBondBond, CrossBondAngle, LinearAngle>;
using SupportedFourCenterTypes =
        TypeList<ProperDihedral, ImproperDihedral, ImproperProperDihedral, RyckaertBellemanDihedral>;
using SupportedFiveCenterTypes = TypeList<Default5Center>;


//! \brief listed types with charges
using PolarizationTypes = TypeList<SimplePolarization>;

/*! \brief pairs with Van-der-Waals and Coulomb interactions
 *
 * Pairs have to be distinguished from other interactions with charges, because they describe
 * non-bonded forces whose energies are accounted separately from the other listed energies.
 */
using PairListedTypes = TypeList<PairLJType, PairLJChargeType>;

//! \brief restraint types
using RestraintTypes = TypeList<PositionRestraints>;

using AggregateTypes =
        TypeList<FourCenterAggregate<HarmonicBondType, HarmonicAngle, PairLJType, ProperDihedral>,
                 FourCenterAggregate<HarmonicBondType, HarmonicAngle, PairLJType, RyckaertBellemanDihedral>,
                 FourCenterAggregate<HarmonicBondType, HarmonicAngle, PairLJType, ImproperDihedral>>;

//***********************************************************************************

using BasicListedTypes =
        Fuse<SupportedTwoCenterTypes, SupportedThreeCenterTypes, SupportedFourCenterTypes, SupportedFiveCenterTypes>;

using ChargedListedTypes = Fuse<PolarizationTypes, PairListedTypes, AggregateTypes>;

using GpuListedTypes = Fuse<BasicListedTypes, PairListedTypes, AggregateTypes>;

using AllListedTypes = Fuse<BasicListedTypes, ChargedListedTypes, RestraintTypes>;

//! \brief meta function to map from an Interaction type to the number of interaction centers
template<class Interaction, class = void>
struct NCenter
{
};

//! \brief meta function return value for restraint interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, RestraintTypes>{}>> :
    util::integral_constant<std::size_t, 1>
{
};

//! \brief meta function return value for two-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedTwoCenterTypes>{}>> :
    util::integral_constant<std::size_t, 2>
{
};

template<>
struct NCenter<SimplePolarization> : util::integral_constant<std::size_t, 2>
{
};

template<>
struct NCenter<PairLJType> : util::integral_constant<std::size_t, 2>
{
};

template<>
struct NCenter<PairLJChargeType> : util::integral_constant<std::size_t, 2>
{
};

//! \brief meta function return value for three-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedThreeCenterTypes>{}>> :
    util::integral_constant<std::size_t, 3>
{
};

template<class Bond, class Angle>
struct NCenter<ThreeCenterAggregate<Bond, Angle>> : util::integral_constant<std::size_t, 3>
{
};

//! \brief meta function return value for four-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedFourCenterTypes>{}>> :
    util::integral_constant<std::size_t, 4>
{
};

template<class... Ts>
struct NCenter<FourCenterAggregate<Ts...>> : util::integral_constant<std::size_t, 4>
{
};

//! \brief meta function return value for five-center interactions
template<class Interaction>
struct NCenter<Interaction, std::enable_if_t<Contains<Interaction, SupportedFiveCenterTypes>{}>> :
    util::integral_constant<std::size_t, 5>
{
};

template<size_t N>
using IndexArray = util::array<int, N>;

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

    // vector of unique interaction parameters for states A & B
    std::vector<InteractionType> parametersA;
    std::vector<InteractionType> parametersB;
    // tuple format: <particleID i, particleID j, ..., InteractionInstanceIndex>
    std::vector<InteractionIndex<InteractionType>> indices;
};

using TwoCenterInteraction   = Reduce<std::variant, SupportedTwoCenterTypes>;
using ThreeCenterInteraction = Reduce<std::variant, SupportedThreeCenterTypes>;
using FourCenterInteraction  = Reduce<std::variant, SupportedFourCenterTypes>;
using FiveCenterInteraction  = Reduce<std::variant, SupportedFiveCenterTypes>;

//! This is the complete type that holds all listed interaction data
// result: std::tuple<ListedTypeData<SupportedListedType1>, ...>
using ListedInteractionData = Reduce<std::tuple, Map<ListedTypeData, AllListedTypes>>;

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_DEFINITIONS_H
