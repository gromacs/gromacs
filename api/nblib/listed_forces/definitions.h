/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \inpublicapi \file
 * \brief
 * Definitions for supported nblib listed interaction data, such as bonds, angles, dihedrals, etc
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * A note on the preprocessor (PP) usage in this file:
 *
 * The PP macros defined here are used exclusively to generate
 * template instantiations declarations of the form "extern template function(X)"
 * in header files and "template function(X)" in .cpp files.
 * These declarations do not affect the program logic in any way and neither are they
 * required to read and understand the behavior of the code as they do not
 * result in any executable instructions.
 * In fact, it would even be technically possible to omit these PP generated
 * declarations in the header files and replace them with an unused static function
 * in the .cpp file that calls the template function in question
 * (e.g. Molecule::addInteraction) once with each type from the variadic template
 * TypeLists declared in this file. This would be enough to create the required instantiations.
 * It would, however, create more work for the compiler which then has to instantiate the
 * templates in the header in each translation unit where the header is included.
 * Doing this results in a compiler warning.
 *
 */
#ifndef NBLIB_LISTEDFORCES_DEFINITIONS_H
#define NBLIB_LISTEDFORCES_DEFINITIONS_H

#include "nblib/util/user.h"
#include "bondtypes.h"

namespace nblib
{

//***********************************************************************************

/*! \brief These macros define what interaction types are supported in
 *  -Molecule
 *  -Topology
 *  -ListedForceCalculator
 *
 *  To enable force calculation for your new interaction type that you've added to bondtypes.h,
 *  list your new type here under the appropriate category and make sure that you've added
 *  a kernel in kernels.hpp
 */

#define SUPPORTED_TWO_CENTER_TYPES \
    HarmonicBondType, G96BondType, CubicBondType, FENEBondType, HalfAttractiveQuarticBondType

#define SUPPORTED_THREE_CENTER_TYPES DefaultAngle

#define SUPPORTED_FOUR_CENTER_TYPES ProperDihedral, ImproperDihedral, RyckaertBellemanDihedral

#define SUPPORTED_FIVE_CENTER_TYPES Default5Center

//***********************************************************************************

#define SUPPORTED_LISTED_TYPES                                                             \
    SUPPORTED_TWO_CENTER_TYPES, SUPPORTED_THREE_CENTER_TYPES, SUPPORTED_FOUR_CENTER_TYPES, \
            SUPPORTED_FIVE_CENTER_TYPES

#define NBLIB_ALWAYS_INLINE __attribute((always_inline))

//! \brief encodes the number of integers needed to represent 2-center interactions (bonds, pairs)
using TwoCenterInteractionIndex = std::array<int, 3>;
//! \brief encodes the number of integers needed to represent 3-center interactions (angles)
using ThreeCenterInteractionIndex = std::array<int, 4>;
//! \brief encodes the number of integers needed to represent 4-center interactions (dihedrals)
using FourCenterInteractionIndex = std::array<int, 5>;
//! \brief encodes the number of integers needed to represent 5-center interactions (CMAP)
using FiveCenterInteractionIndex = std::array<int, 6>;

//! \brief data type for pairwise interactions, e.g. bonds
template<class TwoCenterType>
struct TwoCenterData
{
    using type = TwoCenterType;

    // tuple format: <particleID i, particleID j, TwoCenterInstanceIndex>
    std::vector<TwoCenterInteractionIndex> indices;
    // vector of unique TwoCenterType instances
    std::vector<TwoCenterType> parameters;
};

//! \brief data type for three-center interactions, e.g. angles
template<class ThreeCenterType>
struct ThreeCenterData
{
    using type = ThreeCenterType;

    // tuple format: <particleID i, particleID j, particleID k, ThreeCenterInstanceIndex>
    std::vector<ThreeCenterInteractionIndex> indices;
    // vector of unique ThreeCenterType instances
    std::vector<ThreeCenterType> parameters;
};

//! \brief data type for four-center interactions, e.g. dihedrals
template<class FourCenterType>
struct FourCenterData
{
    using type = FourCenterType;

    // tuple format: <particleID i, particleID j, particleID k, particleID l, FourCenterInstanceIndex>
    std::vector<FourCenterInteractionIndex> indices;
    // vector of unique FiveCenterType instances
    std::vector<FourCenterType> parameters;
};

//! \brief data type for five-center interactions, e.g. CMAP
template<class FiveCenterType>
struct FiveCenterData
{
    using type = FiveCenterType;

    // tuple format: <particleID i, particleID j, particleID k, particleID l, particleID m, FiveCenterInstanceIndex>
    std::vector<FiveCenterInteractionIndex> indices;
    // vector of unique FiveCenterType instances
    std::vector<FiveCenterType> parameters;
};


using SupportedTwoCenterTypes = TypeList<SUPPORTED_TWO_CENTER_TYPES>;
// std::tuple<TwoCenterData<TwoCenterType1>, ...>
using TwoCenterInteractionData = Reduce<std::tuple, Map<TwoCenterData, SupportedTwoCenterTypes>>;

using SupportedThreeCenterTypes = TypeList<SUPPORTED_THREE_CENTER_TYPES>;
// std::tuple<AngleData<ThreeCenterType1>, ...>
using ThreeCenterInteractionData = Reduce<std::tuple, Map<ThreeCenterData, SupportedThreeCenterTypes>>;

using SupportedFourCenterTypes = TypeList<SUPPORTED_FOUR_CENTER_TYPES>;
// std::tuple<FourCenterData<FourCenterType1>, ...>
using FourCenterInteractionData = Reduce<std::tuple, Map<FourCenterData, SupportedFourCenterTypes>>;

using SupportedFiveCenterTypes = TypeList<SUPPORTED_FIVE_CENTER_TYPES>;
// std::tuple<FiveCenterData<FiveCenterType1>, ...>
using FiveCenterInteractionData = Reduce<std::tuple, Map<FiveCenterData, SupportedFiveCenterTypes>>;

//! This is the complete type that holds all listed interaction data
using ListedInteractionData = decltype(std::tuple_cat(TwoCenterInteractionData{},
                                                      ThreeCenterInteractionData{},
                                                      FourCenterInteractionData{},
                                                      FiveCenterInteractionData{}));
} // namespace nblib
#endif // NBLIB_LISTEDFORCES_DEFINITIONS_H
