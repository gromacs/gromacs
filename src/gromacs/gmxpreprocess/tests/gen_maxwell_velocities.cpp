/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Tests for velocity generation.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/gen_maxwell_velocities.h"

#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"
#include "testutils/topologyhelpers.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief
 * Test params for testing velocity generation.
 *
 * Order is: temperature, seed, numWaters
 */
using MaxwellTestParams = std::tuple<real, int, int>;

class MaxwellTest : public ::testing::Test, public ::testing::WithParamInterface<MaxwellTestParams>
{

public:
    MaxwellTest() : checker_(data_.rootChecker()) {}

    ~MaxwellTest() override;
    //! Initialize topology with \p numWaters.
    void initMtop(int numWaters);
    //! Run test with specific \p temp and \p seed.
    void runTest(real temp, int seed);

private:
    //! System topology.
    gmx_mtop_t mtop_;
    //! Velocity vector.
    std::vector<gmx::RVec> v;
    //! Storage for reference data.
    TestReferenceData data_;
    //! Checker for reference data.
    TestReferenceChecker checker_;
};

MaxwellTest::~MaxwellTest()
{
    done_atom(&mtop_.moltype[0].atoms);
}

void MaxwellTest::initMtop(int numWaters)
{
    addNWaterMolecules(&mtop_, numWaters);
    mtop_.finalize();
    v.resize(mtop_.natoms);
}

void MaxwellTest::runTest(real temp, int seed)
{
    MDLogger logger;
    maxwell_speed(temp, seed, &mtop_, as_rvec_array(v.data()), logger);
    TestReferenceChecker compound(checker_.checkCompound("Velocities", nullptr));
    const auto           tolerance = relativeToleranceAsPrecisionDependentUlp(1.0, 40, 20);
    compound.setDefaultTolerance(tolerance);
    compound.checkSequence(v.begin(), v.end(), "Velocity values");
}

TEST_P(MaxwellTest, CreationWorks)
{
    const auto& params    = GetParam();
    const real  temp      = std::get<0>(params);
    const int   seed      = std::get<1>(params);
    const int   numWaters = std::get<2>(params);
    initMtop(numWaters);

    runTest(temp, seed);
}

INSTANTIATE_TEST_SUITE_P(CorrectVelocity,
                         MaxwellTest,
                         ::testing::Combine(::testing::Values(150, 298, 313, 350),
                                            ::testing::Values(1, 42),
                                            ::testing::Values(23, 42)));

} // namespace
} // namespace test
} // namespace gmx
