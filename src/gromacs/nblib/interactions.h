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
/*! \file
 * \brief
 * Implements nblib particle-types interactions
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#ifndef GMX_NBLIB_INTERACTIONS_H
#define GMX_NBLIB_INTERACTIONS_H

#include <map>
#include <set>

#include "gromacs/math/vectypes.h"

#include "particletype.h"

namespace nblib
{
using ParticleTypeName        = std::string;
using NonBondedInteractionMap = std::map<std::tuple<ParticleTypeName, ParticleTypeName>, std::tuple<C6, C12>>;

//! Enum for selecting coulomb type
enum class CombinationRule : int
{
    Geometric,
    Count
};

namespace detail
{

real combineNonbondedParameters(real v, real w, CombinationRule combinationRule);

}


class ParticleTypesInteractions
{
public:
    ParticleTypesInteractions() = default;

    void add(ParticleType particleType, C6 c6, C12 c12);

    void add(ParticleType particleType1, ParticleType particleType2, C6 c6, C12 c12);

    NonBondedInteractionMap generateTable(CombinationRule combinationRule);

private:

    std::unordered_map<ParticleTypeName, std::tuple<C6, C12>> singleParticleInteractionsMap_;
    std::map<std::tuple<ParticleTypeName, ParticleTypeName>, std::tuple<C6, C12>> twoParticlesInteractionsMap_;

    // Helper data structure to collect the unique ParticleTypes
    std::set<ParticleTypeName> particleTypesSet_;
};

} // namespace nblib
#endif // GMX_NBLIB_INTERACTIONS_H
