/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests structure similarity measures rmsd and size-independent rho factor.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/do_fit.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"

namespace
{

using gmx::RVec;
using gmx::test::defaultRealTolerance;

TEST(StructureSimilarity, StructureComparedToSelfHasZeroRMSD)
{
    std::vector<RVec>       structure = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector<real>       masses    = {1, 1, 1};
    auto                    x1        = gmx::as_rvec_array(structure.data());

    const real              rmsdSameStructure = rmsdev(masses.size(), masses.data(), x1, x1);

    EXPECT_REAL_EQ_TOL(0., rmsdSameStructure, defaultRealTolerance());
}

TEST(StructureSimilarity, StructureComparedToSelfHasZeroRho)
{
    std::vector<RVec>       structure = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector<real>       masses    = {1, 1, 1};
    auto                    x1        = gmx::as_rvec_array(structure.data());

    const real              rmsdSameStructure = rhodev(masses.size(), masses.data(), x1, x1);

    EXPECT_REAL_EQ_TOL(0., rmsdSameStructure, defaultRealTolerance());
}

TEST(StructureSimilarity, YieldsCorrectRMSD)
{
    std::vector<RVec>       structureA = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector<RVec>       structureB = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
    std::vector<real>       masses     = {1, 1, 1};
    auto                    x1         = gmx::as_rvec_array(structureA.data());
    auto                    x2         = gmx::as_rvec_array(structureB.data());

    const real              rmsdSameStructure = rmsdev(masses.size(), masses.data(), x1, x2);
    constexpr real          sqrtTwo           = 1.41421356237;
    EXPECT_REAL_EQ_TOL(sqrtTwo, rmsdSameStructure, defaultRealTolerance());
}

TEST(StructureSimilarity, YieldsCorrectRho)
{
    std::vector<RVec>       structureA = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector<RVec>       structureB = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
    std::vector<real>       masses     = {1, 1, 1};
    const auto              x1         = gmx::as_rvec_array(structureA.data());
    const auto              x2         = gmx::as_rvec_array(structureB.data());

    const real              rhoSameStructure = rhodev(masses.size(), masses.data(), x1, x2);

    EXPECT_REAL_EQ_TOL(2., rhoSameStructure, defaultRealTolerance());
}

TEST(StructureSimilarity, YieldsCorrectRMSDWithIndex)
{
    std::vector<RVec>       structureA = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
    std::vector<RVec>       structureB = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {0, 0, 0}};
    std::vector<real>       masses     = {1, 1, 1, 0};
    std::vector<int>        index      = {0, 1, 2};
    auto                    x1         = gmx::as_rvec_array(structureA.data());
    auto                    x2         = gmx::as_rvec_array(structureB.data());

    const real              rmsdSameStructure = rmsdev_ind(index.size(), index.data(), masses.data(), x1, x2);
    constexpr real          sqrtTwo           = 1.41421356237;
    EXPECT_REAL_EQ_TOL(sqrtTwo, rmsdSameStructure, defaultRealTolerance());
}

TEST(StructureSimilarity, YieldsCorrectRhoWidthIndex)
{
    std::vector<RVec>       structureA = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
    std::vector<RVec>       structureB = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {0, 0, 0}};
    std::vector<real>       masses     = {1, 1, 1, 1};
    std::vector<int>        index      = {0, 1, 2};
    auto                    x1         = gmx::as_rvec_array(structureA.data());
    auto                    x2         = gmx::as_rvec_array(structureB.data());

    const real              rhoSameStructure = rhodev_ind(index.size(), index.data(), masses.data(), x1, x2);

    EXPECT_REAL_EQ_TOL(2., rhoSameStructure, defaultRealTolerance());
}


} // namespace
