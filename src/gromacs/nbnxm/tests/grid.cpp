/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for the NBNxM search grid functionality.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/nbnxm/grid.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/nbnxm/gridsetdata.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/utility/enumerationhelpers.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

constexpr int sc_numThreads = 1;

enum class CoordDistribution : int
{
    Uniform, //!< Fully uniform distribution over the whole volume
    Column,  //!< A column of particles extending over the whole z-range
    Layer,   //!< A layer of particles extending over the whole x/y-range
    Block,   //!< A block of particles a corner
    Count
};

//! Lower corner for the grids
static const RVec sc_lowerCorner = { 0, 0, 0 };

//! Upper corners for the grids
static const EnumerationArray<CoordDistribution, RVec> sc_upperCorner = {
    { { 2, 3, 5 }, { 16, 8, 1 }, { 2, 2, 20 }, { 4, 4, 4 } }
};

//! A list of 4 particle coordinates that fill a grid cell with uniform x, y and z-distribution
static const std::array<RVec, 4> sc_cellX{ { { 0.25_real, 0.25_real, 0.25_real },
                                             { 0.25_real, 0.75_real, 0.25_real },
                                             { 0.25_real, 0.25_real, 0.75_real },
                                             { 0.25_real, 0.75_real, 0.75_real } } };

//! Returns coordinates for the given setup
std::vector<RVec> generateCoords(const CoordDistribution coordDistribution)
{
    // Note that to keep the atom count low here we need to choose the range
    // where we put particles such that the boundary coincides with the cell
    // boundaries of the initially chosen search grid.
    BasicVector<Range<int>> range;

    switch (coordDistribution)
    {
        case CoordDistribution::Uniform:
            for (int d = 0; d < DIM; d++)
            {
                range[d] = { 0, int(sc_upperCorner[coordDistribution][d]) };
            }
            break;
        case CoordDistribution::Column:
            range[XX] = { 4, 8 };
            range[YY] = { 4, 8 };
            range[ZZ] = { 0, int(sc_upperCorner[coordDistribution][ZZ]) };
            break;
        case CoordDistribution::Layer:
            for (int d = 0; d < ZZ; d++)
            {
                range[d] = { 0, int(sc_upperCorner[coordDistribution][d]) };
            }
            range[ZZ] = { 2, 4 };
            break;
        case CoordDistribution::Block:
            range[XX] = { 2, 4 };
            range[YY] = { 2, 4 };
            range[ZZ] = { 0, 2 };
            break;
        case CoordDistribution::Count: GMX_RELEASE_ASSERT(false, "We should not use Count"); break;
    }

    std::vector<RVec> x;

    for (int i : range[XX])
    {
        for (int j : range[YY])
        {
            for (int k : range[ZZ])
            {
                for (const RVec& cx : sc_cellX)
                {
                    x.push_back({ i + cx[XX], j + cx[YY], k + cx[ZZ] });
                }
            }
        }
    }

    return x;
}

//! Class that sets up and holds a search grid
class GridTest : public ::testing::TestWithParam<CoordDistribution>
{
    static constexpr bool sc_haveFep = false;

public:
    GridTest() :
        grid_(PairlistType::Simple4x4, 0, sc_haveFep, HostAllocationPolicy{}),
        gridWork_(sc_numThreads),
        x_(generateCoords(GetParam()))
    {
        gmx_omp_nthreads_set(ModuleMultiThread::Pairsearch, sc_numThreads);
    }

    //! Generates the search grid, return the computed uniform and effective atom densities
    std::pair<real, real> generateGrid()
    {
        const int numAtoms = gmx::ssize(x_);

        real uniformAtomDensity = -1;

        const real effectiveAtomDensity = grid_.generateAndFill2D(gridWork_,
                                                                  &bins_,
                                                                  sc_lowerCorner,
                                                                  sc_upperCorner[GetParam()],
                                                                  nullptr,
                                                                  { 0, numAtoms },
                                                                  numAtoms,
                                                                  &uniformAtomDensity,
                                                                  0,
                                                                  x_,
                                                                  0,
                                                                  nullptr,
                                                                  true);

        return { uniformAtomDensity, effectiveAtomDensity };
    }

    //! Returns the uniform atom density
    real getDensity() const
    {
        return gmx::ssize(x_)
               / (sc_upperCorner[GetParam()][XX] * sc_upperCorner[GetParam()][YY]
                  * sc_upperCorner[GetParam()][ZZ]);
    }

private:
    Grid                  grid_;
    std::vector<GridWork> gridWork_;
    HostVector<int>       bins_;
    std::vector<RVec>     x_;
};

TEST_P(GridTest, CheckDensities)
{
    // The density does not need to be accurate. We want to avoid being off by a factor of >>2.
    // With the low number of particles in this test we get a bit below 15%.
    const FloatingPointTolerance tol = relativeToleranceAsFloatingPoint(1, 0.15);

    real uniformAtomDensity;
    real effectiveAtomDensity;
    std::tie(uniformAtomDensity, effectiveAtomDensity) = generateGrid();

    const real refUniformAtomDensity = getDensity();

    EXPECT_REAL_EQ_TOL(refUniformAtomDensity, uniformAtomDensity, tol);
    EXPECT_REAL_EQ_TOL(sc_cellX.size(), effectiveAtomDensity, tol);
}

INSTANTIATE_TEST_SUITE_P(WithParameters,
                         GridTest,
                         ::testing::Values(CoordDistribution::Uniform,
                                           CoordDistribution::Column,
                                           CoordDistribution::Layer,
                                           CoordDistribution::Block));

} // namespace
} // namespace test
} // namespace gmx
