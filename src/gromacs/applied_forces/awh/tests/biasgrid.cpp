/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/applied_forces/awh/biasgrid.h"

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/applied_forces/awh/dimparams.h"
#include "gromacs/applied_forces/awh/tests/awh_setup.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

TEST(biasGridTest, neighborhood)
{
    constexpr double pointsPerScope = BiasGrid::c_scopeCutoff * BiasGrid::c_numPointsPerSigma;
    GMX_RELEASE_ASSERT(std::abs(pointsPerScope - std::round(pointsPerScope)) > 1e-4,
                       "If the scope is close to an integer number of points, this test can be "
                       "unstable due to rounding issues");

    const int scopeInPoints = static_cast<int>(pointsPerScope);

    /* We test a 2-dimensional grid with dim0 periodic, dim1 not periodic.
     * This test should cover most relevant cases: multi-dimension grids
     * and non-periodic and periodic dimensions.
     * The only thing not tested here is a periodic dimension with a gap.
     */

    const int                 numDim = 2;
    std::vector<AwhDimParams> awhDimParams;
    AwhCoordinateProviderType coordinateProvider = AwhCoordinateProviderType::Pull;
    double                    diffusion          = 0.1;
    {
        int    coordIndex = 0;
        double origin     = -5;
        double end        = 5;
        double period     = 10;
        auto   awhDimBuffer =
                awhDimParamSerialized(coordinateProvider, coordIndex, origin, end, period, diffusion);
        gmx::InMemoryDeserializer serializer(awhDimBuffer, false);
        awhDimParams.emplace_back(&serializer);
    }
    {
        int    coordIndex = 1;
        double origin     = 0.5;
        double end        = 2;
        double period     = 0;
        auto   awhDimBuffer =
                awhDimParamSerialized(coordinateProvider, coordIndex, origin, end, period, diffusion);
        gmx::InMemoryDeserializer serializer(awhDimBuffer, false);
        awhDimParams.emplace_back(&serializer);
    }

    const real conversionFactor = 1;
    const real beta             = 3.0;

    /* Set up dimParams to get about 15 points along each dimension */
    std::vector<DimParams> dimParams;
    dimParams.push_back(DimParams::pullDimParams(conversionFactor, 1 / (beta * 0.7 * 0.7), beta));
    dimParams.push_back(DimParams::pullDimParams(conversionFactor, 1 / (beta * 0.1 * 0.1), beta));

    BiasGrid grid(dimParams, awhDimParams);

    const int numPoints = grid.numPoints();

    int numPointsDim[numDim];
    for (int d = 0; d < numDim; d++)
    {
        numPointsDim[d] = grid.axis(d).numPoints();

        GMX_RELEASE_ASSERT(numPointsDim[d] >= 2 * scopeInPoints + 1,
                           "The test code currently assume that the grid is larger of equal to the "
                           "scope in each dimension");
    }
    int halfNumPoints0 = numPointsDim[0] / 2;

    bool haveOutOfGridNeighbors  = false;
    bool haveDuplicateNeighbors  = false;
    bool haveIncorrectNeighbors  = false;
    bool haveCorrectNumNeighbors = true;

    /* Set up a grid for checking for duplicate neighbors */
    std::vector<bool> isInNeighborhood(grid.numPoints(), false);

    /* Checking for all points is overkill, we check every 7th */
    for (size_t i = 0; i < grid.numPoints(); i += 7)
    {
        const GridPoint& point = grid.point(i);

        /* NOTE: This code relies on major-minor index ordering in Grid */
        int pointIndex0 = i / numPointsDim[1];
        int pointIndex1 = i - pointIndex0 * numPointsDim[1];

        /* To check if we have the correct neighbors, we check the expected
         * number of neighbors, if all neighbors are within the grid bounds
         * and if they are within scope.
         */
        int    distanceFromEdge1 = std::min(pointIndex1, numPointsDim[1] - 1 - pointIndex1);
        size_t numNeighbors      = (2 * scopeInPoints + 1)
                              * (scopeInPoints + std::min(scopeInPoints, distanceFromEdge1) + 1);
        if (point.neighbor.size() != numNeighbors)
        {
            haveCorrectNumNeighbors = false;
        }

        for (const auto& j : point.neighbor)
        {
            if (j >= 0 && j < numPoints)
            {
                if (isInNeighborhood[j])
                {
                    haveDuplicateNeighbors = true;
                }
                isInNeighborhood[j] = true;

                int neighborIndex0 = j / numPointsDim[1];
                int neighborIndex1 = j - neighborIndex0 * numPointsDim[1];

                int distance0 = neighborIndex0 - pointIndex0;
                int distance1 = neighborIndex1 - pointIndex1;
                /* Adjust distance for periodicity of dimension 0 */
                if (distance0 < -halfNumPoints0)
                {
                    distance0 += numPointsDim[0];
                }
                else if (distance0 > halfNumPoints0)
                {
                    distance0 -= numPointsDim[0];
                }
                /* Check if the distance is within scope */
                if (distance0 < -scopeInPoints || distance0 > scopeInPoints
                    || distance1 < -scopeInPoints || distance1 > scopeInPoints)
                {
                    haveIncorrectNeighbors = true;
                }
            }
            else
            {
                haveOutOfGridNeighbors = true;
            }
        }

        /* Clear the marked points in the checking grid */
        for (const auto& neighbor : point.neighbor)
        {
            if (neighbor >= 0 && neighbor < numPoints)
            {
                isInNeighborhood[neighbor] = false;
            }
        }
    }

    EXPECT_FALSE(haveOutOfGridNeighbors);
    EXPECT_FALSE(haveDuplicateNeighbors);
    EXPECT_FALSE(haveIncorrectNeighbors);
    EXPECT_TRUE(haveCorrectNumNeighbors);
}

} // namespace test
} // namespace gmx
