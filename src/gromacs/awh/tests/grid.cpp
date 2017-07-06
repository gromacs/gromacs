/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/awh/grid.h"

#include <cmath>

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/utility/gmxassert.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

TEST(gridTest, neighborhood)
{
    constexpr double pointsPerScope = Grid::c_scopeCutoff*Grid::c_numPointsPerSigma;
    GMX_RELEASE_ASSERT(std::abs(pointsPerScope - std::round(pointsPerScope)) > 1e-4, "If the scope is close to an integer number of points, this test can be unstable due to rounding issues");

    const int scopeInPoints = static_cast<int>(pointsPerScope);

    /* We test a 2-dimensional grid with dim0 periodic, dim1 not periodic.
     * This test should cover most relevant cases: multi-dimension grids
     * and non-periodic and periodic dimensions.
     * The only thing not tested here is a periodic dimension with a gap.
     */

    const int                 numDim = 2;
    std::vector<AwhDimParams> awhDimParams(numDim);

    awhDimParams[0].origin      = -5;
    awhDimParams[0].end         = 5;
    awhDimParams[0].period      = 10;

    awhDimParams[1].origin      = 0.5;
    awhDimParams[1].end         = 2.0;
    awhDimParams[1].period      =   0;

    const real conversionFactor = 1;
    const real beta             = 3.0;

    /* Set up dimParams to get about 15 points along each dimension */
    std::vector<DimParams> dimParams;
    dimParams.push_back(DimParams(conversionFactor, 1/(beta*0.7*0.7), beta));
    dimParams.push_back(DimParams(conversionFactor, 1/(beta*0.1*0.1), beta));

    Grid       grid(dimParams, awhDimParams.data());

    const int  numPoints         = grid.numPoints();

    int        numPointsDim[numDim];
    for (int d = 0; d < numDim; d++)
    {
        numPointsDim[d]          = grid.axis(d).numPoints();

        GMX_RELEASE_ASSERT(numPointsDim[d] >= 2*scopeInPoints + 1, "The test code currently assume that the grid is larger of equal to the scope in each dimension");
    }
    int  halfNumPoints0          = numPointsDim[0]/2;

    bool haveOutOfGridNeighbors  = false;
    bool haveDuplicateNeighbors  = false;
    bool haveIncorrectNeighbors  = false;
    bool haveCorrectNumNeighbors = true;

    /* Set up a grid for checking for duplicate neighbors */
    std::vector<bool> isInNeighborhood(grid.numPoints(), false);

    /* Checking for all points is overkill, we check every 7th */
    for (size_t i = 0; i < grid.numPoints(); i += 7)
    {
        const GridPoint &point = grid.point(i);

        /* NOTE: This code relies on major-minor index ordering in Grid */
        int    pointIndex0       = i/numPointsDim[1];
        int    pointIndex1       = i - pointIndex0*numPointsDim[1];

        /* To check if we have the correct neighbors, we check the expected
         * number of neighbors, if all neighbors are within the grid bounds
         * and if they are within scope.
         */
        int    distanceFromEdge1 = std::min(pointIndex1, numPointsDim[1] - 1 - pointIndex1);
        size_t numNeighbors      = (2*scopeInPoints + 1)*(scopeInPoints + std::min(scopeInPoints, distanceFromEdge1) + 1);
        if (point.neighbor.size() != numNeighbors)
        {
            haveCorrectNumNeighbors = false;
        }

        for (auto &j : point.neighbor)
        {
            if (j >= 0 && j < numPoints)
            {
                if (isInNeighborhood[j])
                {
                    haveDuplicateNeighbors = true;
                }
                isInNeighborhood[j] = true;

                int neighborIndex0  = j/numPointsDim[1];
                int neighborIndex1  = j - neighborIndex0*numPointsDim[1];

                int distance0       = neighborIndex0 - pointIndex0;
                int distance1       = neighborIndex1 - pointIndex1;
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
                if (distance0 < -scopeInPoints || distance0 > scopeInPoints ||
                    distance1 < -scopeInPoints || distance1 > scopeInPoints)
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
        for (auto &neighbor : point.neighbor)
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

} // namespace
} // namespace
