/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests structure similarity measures rmsd and size-independent rho factor.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include <cmath>

#include <array>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::RVec;
using gmx::test::defaultRealTolerance;
class StructureSimilarityTest : public ::testing::Test
{
protected:
    static constexpr int       c_nAtoms = 4;
    std::array<RVec, c_nAtoms> structureA_{ { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 } } };
    std::array<RVec, c_nAtoms> structureB_{ { { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 0 }, { 0, 0, 0 } } };
    std::array<real, c_nAtoms> masses_{ { 1, 1, 1, 0 } };
    std::array<int, 3>         index_{ { 0, 1, 2 } };
    rvec*                      x1_ = gmx::as_rvec_array(structureA_.data());
    rvec*                      x2_ = gmx::as_rvec_array(structureB_.data());
    real*                      m_  = masses_.data();
};

TEST_F(StructureSimilarityTest, StructureComparedToSelfHasZeroRMSD)
{
    EXPECT_REAL_EQ_TOL(0., rmsdev(c_nAtoms, m_, x1_, x1_), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, StructureComparedToSelfHasZeroRho)
{
    EXPECT_REAL_EQ_TOL(0., rhodev(c_nAtoms, m_, x1_, x1_), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRMSD)
{
    EXPECT_REAL_EQ_TOL(std::sqrt(2.0), rmsdev(c_nAtoms, m_, x1_, x2_), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRho)
{
    EXPECT_REAL_EQ_TOL(2., rhodev(c_nAtoms, m_, x1_, x2_), defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRMSDWithIndex)
{
    EXPECT_REAL_EQ_TOL(std::sqrt(2.0),
                       rmsdev_ind(index_.size(), index_.data(), m_, x1_, x2_),
                       defaultRealTolerance());
}

TEST_F(StructureSimilarityTest, YieldsCorrectRhoWidthIndex)
{
    EXPECT_REAL_EQ_TOL(2., rhodev_ind(index_.size(), index_.data(), m_, x1_, x2_), defaultRealTolerance());
}

} // namespace
} // namespace test
} // namespace gmx
