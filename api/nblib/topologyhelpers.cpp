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
 * Implements nblib Topology helpers
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#include "topologyhelpers.h"

#include <algorithm>

#include "gromacs/topology/exclusionblocks.h"

namespace nblib
{

std::vector<gmx::ExclusionBlock> toGmxExclusionBlock(const std::vector<std::tuple<int, int>>& tupleList)
{
    std::vector<gmx::ExclusionBlock> ret;

    auto firstLowerThan = [](auto const& tup1, auto const& tup2) {
        return std::get<0>(tup1) < std::get<0>(tup2);
    };

    // Note this is a programming error as all particles should exclude at least themselves and empty topologies are not allowed.
    // Note also that this is also checked in the parent function with a more informative error message.
    assert((!tupleList.empty() && "No exclusions found.\n"));

    // initialize pair of iterators delimiting the range of exclusions for
    // the first particle in the list
    auto range = std::equal_range(std::begin(tupleList), std::end(tupleList), tupleList[0], firstLowerThan);
    auto it1 = range.first;
    auto it2 = range.second;

    // loop over all exclusions in molecule, linear in tupleList.size()
    while (it1 != std::end(tupleList))
    {
        gmx::ExclusionBlock localBlock;
        // loop over all exclusions for current particle
        for (; it1 != it2; ++it1)
        {
            localBlock.atomNumber.push_back(std::get<1>(*it1));
        }

        ret.push_back(localBlock);

        // update the upper bound of the range for the next particle
        if (it1 != end(tupleList))
        {
            it2 = std::upper_bound(it1, std::end(tupleList), *it1, firstLowerThan);
        }
    }

    return ret;
}

std::vector<gmx::ExclusionBlock> offsetGmxBlock(std::vector<gmx::ExclusionBlock> inBlock, int offset)
{
    // shift particle numbers by offset
    for (auto& localBlock : inBlock)
    {
        std::transform(std::begin(localBlock.atomNumber),
                       std::end(localBlock.atomNumber),
                       std::begin(localBlock.atomNumber),
                       [offset](auto i) { return i + offset; });
    }

    return inBlock;
}

} // namespace nblib
