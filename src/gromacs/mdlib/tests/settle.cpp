/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
 * \brief Tests for SETTLE constraints
 *
 * The test runs on several small systems, containing 1 to 17 water molecules,
 * with and without periodic boundary conditions, with and without velocity
 * and virial updates. The CPU and GPU versions are tested, if the code was
 * compiled with CUDA support and there is a CUDA-capable GPU in the system.
 *
 * The tests check:
 * 1. If the final distances between constrained atoms are within tolerance
 *    from the target distance.
 * 2. If the velocities were updated when needed.
 * 3. If the virial was computed.
 *
 * The test also compares the results from the CPU and GPU versions of the
 * algorithm: final coordinates, velocities and virial should be within
 * tolerance to one another.
 *
 * \todo This also tests that if the calling code requires velocities
 *       and virial updates, that those outputs do change, but does not
 *       test that those changes are correct.
 *
 * \todo Only no-PBC and cubic-PBC are tested here, but the correct
 *       function of the SIMD version of set_pbx_auic in all cases
 *       should be tested elsewhere.
 *
 * \todo The CPU and GPU versions are tested against each other. This
 *       should be changed to a proper test against pre-computed
 *       reference values. Also, these test will dry-run on a CUDA
 *       build if no CUDA-capable GPU is available.
 *
 * \author Mark Abraham  <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/settle.h"

#include "config.h"

#include <algorithm>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/gpu_testutils.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "gromacs/mdlib/tests/watersystem.h"
#include "testutils/testasserts.h"

#include "settletestdata.h"
#include "settletestrunners.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Parameters that will vary from test to test.
 *
 * The test will run in the following parameter space:
 * 1. Number of water molecules (SETTLEs) [1, 2, 4, 5, 7, 10, 12, 15, 17]
 * 2. If the velocities should be updated while constraining [true/false]
 * 3. If the virial should be computed [true/false]
 * 4. Periodic boundary conditions [XYZ/None]
 * 5. Subroutine to use [CPU/GPU version of SETTLE]
 */
typedef std::tuple<int, bool, bool, std::string, std::string> SettleTestParameters;


//! List of all availible algortihms/implementations.
std::vector<std::string> algorithmsNames;

//! Method that fills and returns algorithmNames to the test macros.
std::vector<std::string> getAlgorithmsNames()
{
    algorithmsNames.emplace_back("SETTLE");
    // TODO: Here we should check that at least 1 suitable GPU is available
    if (GMX_GPU == GMX_GPU_CUDA && canComputeOnGpu())
    {
        algorithmsNames.emplace_back("SETTLE_GPU");
    }
    return algorithmsNames;
}

/*! \brief Test fixture for testing SETTLE position updates.
 *
 */
class SettleTest : public ::testing::TestWithParam<SettleTestParameters>
{
    public:
        //! PBC setups
        std::unordered_map <std::string, t_pbc>                     pbcs_;
        //! Algorithms (CPU and GPU versions of SETTLE)
        std::unordered_map <std::string, void(*)(SettleTestData    *testData,
                                                 const t_pbc        pbc,
                                                 const bool         updateVelocities,
                                                 const bool         calcVirial,
                                                 const std::string &testDescription)> algorithms_;

        /*! \brief Test setup function.
         *
         * Setting up the PBCs and algorithms. Note, that corresponding string keywords
         * have to be explicitly added at the end of this file when the tests are called.
         *
         */
        void SetUp() override
        {

            //
            // PBC initialization
            //
            t_pbc pbc;

            // Infinitely small box
            matrix boxNone = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
            set_pbc(&pbc, epbcNONE, boxNone);
            pbcs_["PBCNone"] = pbc;

            // Rectangular box
            matrix boxXyz = {{real(1.86206), 0, 0}, {0, real(1.86206), 0}, {0, 0, real(1.86206)}};
            set_pbc(&pbc, epbcXYZ, boxXyz);
            pbcs_["PBCXYZ"] = pbc;

            //
            // Algorithms
            //
            algorithms_["SETTLE"]     = applySettle;
            algorithms_["SETTLE_GPU"] = applySettleGpu;
        }

        /*! \brief Check if the positions and velocities were constrained correctly.
         *
         * The final distance between constrained atoms is compared to the target distance. The velocities are
         * not properly tested, only the fact that they were updated is checked.
         *
         * \todo Make a proper test for both coordinates and velocities. They should be tested against known
         *       the final values.
         *
         * \param[in]  numSettles        Number of water molecules in the tested system.
         * \param[in]  updateVelocities  If the velocities were updated.
         * \param[in]  tolerance         Tolerance to compare floating point numbers.
         * \param[in]  testDescription   Brief description that will be printed in case of test failure.
         * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
         */
        void checkFinalPositionsAndVelocities(const int                     numSettles,
                                              const bool                    updateVelocities,
                                              const FloatingPointTolerance  tolerance,
                                              const std::string            &testDescription,
                                              const SettleTestData         &testData)
        {
            for (int i = 0; i < numSettles; ++i)
            {
                const gmx::RVec &positionO  = testData.xPrime_[i*3 + 0];
                const gmx::RVec &positionH1 = testData.xPrime_[i*3 + 1];
                const gmx::RVec &positionH2 = testData.xPrime_[i*3 + 2];

                real             dOH = testData.dOH_;
                real             dHH = testData.dHH_;

                EXPECT_REAL_EQ_TOL(dOH*dOH, distance2(positionO, positionH1), tolerance) << formatString("for water %d ", i) << testDescription;
                EXPECT_REAL_EQ_TOL(dOH*dOH, distance2(positionO, positionH2), tolerance) << formatString("for water %d ", i) << testDescription;
                EXPECT_REAL_EQ_TOL(dHH*dHH, distance2(positionH1, positionH2), tolerance) << formatString("for water %d ", i) << testDescription;

                // This merely tests whether the velocities were
                // updated from the starting values of zero (or not),
                // but not whether the update was correct.
                for (int j = 0; j < testData.atomsPerSettle_; ++j)
                {
                    for (int d = 0; d < DIM; ++d)
                    {
                        EXPECT_TRUE(updateVelocities == (0. != testData.v_[i*3 + j][d]))
                        << formatString("for water %d velocity atom %d dim %d", i, j, d)
                        << testDescription;
                    }
                }
            }
        }

        /*! \brief Check if the virial was updated.
         *
         * The two tests on virial are:
         * 1. If it was updated in case calcVirial is true.
         * 2. If it is symmetrical.
         *
         * \todo Test against pre-computed values should be added.
         *
         * \param[in]  calcVirial        If the virial is computed.
         * \param[in]  tolerance         Tolerance to compare floating point numbers.
         * \param[in]  testDescription   Brief description that will be printed in case of test failure.
         * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
         */
        void checkVirial(const bool                   calcVirial,
                         const FloatingPointTolerance tolerance,
                         const std::string           &testDescription,
                         const SettleTestData        &testData)
        {
            // This merely tests whether the viral was updated from
            // the starting values of zero (or not), but not whether
            // any update was correct. The symmetry of the tensor is
            // also tested.
            for (int d = 0; d < DIM; ++d)
            {
                for (int dd = 0; dd < DIM; ++dd)
                {

                    EXPECT_TRUE(calcVirial == (0. != testData.virial_[d][dd]))
                    << formatString("for virial component[%d][%d] ", d, dd)
                    << testDescription;

                    if (calcVirial)
                    {
                        EXPECT_REAL_EQ_TOL(testData.virial_[d][dd], testData.virial_[dd][d], tolerance)
                        << formatString("Virial is not symmetrical for [%d][%d] ", d, dd)
                        << testDescription;
                    }
                }
            }
        }

        /*! \brief Check if two sets of coordinates are equal up to the provided tolerance.
         *
         * This is used to compare the positions and velocities, obtained on the CPU and GPU between
         * each other.
         *
         * \todo Values should be compared to pre-computed reference values instead.
         *
         * \param[in]  coordinates       Coordinates to check.
         * \param[in]  coordinatesRef    Reference values for the constrained coordinates.
         * \param[in]  tolerance         Tolerance to compare floating point numbers
         * \param[in]  testDescription   Brief description that will be printed in case of test failure.
         */
        void checkCoordinatesEqual(const gmx::PaddedVector<gmx::RVec> &coordinates,
                                   const gmx::PaddedVector<gmx::RVec> &coordinatesRef,
                                   const FloatingPointTolerance        tolerance,
                                   const std::string                  &testDescription)
        {
            EXPECT_EQ(coordinates.size(), coordinatesRef.size());
            for (index i = 0; i < ssize(coordinates); ++i)
            {
                for (int d = 0; d < DIM; ++d)
                {
                    EXPECT_REAL_EQ_TOL(coordinates[i][d], coordinatesRef[i][d], tolerance)
                    << formatString("Different coordinates in GPU and CPU codes for particle %zd component %d ", i, d)
                    << testDescription;
                }
            }
        }

        /*! \brief Check if two tensors are equal up to the tolerance.
         *
         * This is used to check if the computation of virial works similarly on CPU and GPU.
         *
         * \todo Values should be compared to pre-computed reference values instead.
         *
         * \param[in]  virial            Tensor to check.
         * \param[in]  virialRef         Reference values for the tensor.
         * \param[in]  tolerance         Tolerance to compare floating point numbers.
         * \param[in]  testDescription   Brief description that will be printed in case of test failure.
         */
        void checkTensorsEqual(const tensor                  virial,
                               const tensor                  virialRef,
                               const FloatingPointTolerance  tolerance,
                               const std::string            &testDescription)
        {
            for (int d1 = 0; d1 < DIM; ++d1)
            {
                for (int d2 = 0; d2 < DIM; ++d2)
                {
                    EXPECT_REAL_EQ_TOL(virial[d1][d2], virialRef[d1][d2], tolerance)
                    << formatString("Different tensor components in GPU and CPU codes for components %d and %d ", d1, d2)
                    << testDescription;
                }
            }
        }

        //! Store whether any compatible GPUs exist.
        static bool s_hasCompatibleGpus;
        //! Before any test is run, work out whether any compatible GPUs exist.
        static void SetUpTestCase()
        {
            s_hasCompatibleGpus = canComputeOnGpu();
        }
};

bool SettleTest::s_hasCompatibleGpus = false;

TEST_P(SettleTest, SatisfiesConstraints)
{
    int         numSettles;
    bool        updateVelocities, calcVirial;
    std::string pbcName;
    std::string algorithmName;
    // Make some symbolic names for the parameter combination.
    std::tie(numSettles, updateVelocities, calcVirial, pbcName, algorithmName) = GetParam();

    // Make a string that describes which parameter combination is
    // being tested, to help make failing tests comprehensible.
    std::string testDescription = formatString("while testing %s with %d SETTLEs, %s, %svelocities and %scalculating the virial",
                                               algorithmName.c_str(),
                                               numSettles,
                                               pbcName.c_str(),
                                               updateVelocities ? "with " : "without ",
                                               calcVirial ? "" : "not ");

    auto testData = std::make_unique<SettleTestData>(numSettles);

    ASSERT_LE(numSettles, testData->xPrime_.size() / testData->atomsPerSettle_) << "cannot test that many SETTLEs " << testDescription;

    t_pbc         pbc = pbcs_.at(pbcName);

    algorithms_.at(algorithmName)(testData.get(),
                                  pbc,
                                  updateVelocities,
                                  calcVirial,
                                  testDescription);

    // The necessary tolerances for the test to pass were determined
    // empirically. This isn't nice, but the required behavior that
    // SETTLE produces constrained coordinates consistent with
    // sensible sampling needs to be tested at a much higher level.
    real                   dOH             = testData->dOH_;
    FloatingPointTolerance tolerance       = relativeToleranceAsPrecisionDependentUlp(dOH*dOH, 80, 380);
    FloatingPointTolerance toleranceVirial = absoluteTolerance(0.000001);

    checkFinalPositionsAndVelocities(numSettles, updateVelocities, tolerance, testDescription, *testData);
    checkVirial(calcVirial, toleranceVirial, testDescription, *testData);
}

#if GMX_GPU == GMX_GPU_CUDA

TEST_P(SettleTest, CompareCpuAndGpu)
{
    // CUDA version will be tested only if:
    // 1. The code was compiled with CUDA
    // 2. There is a CUDA-capable GPU in a system
    // 3. This GPU is detectable
    // 4. GPU detection was not disabled by GMX_DISABLE_GPU_DETECTION environment variable
    if (s_hasCompatibleGpus)
    {
        int         numSettles;
        bool        updateVelocities, calcVirial;
        std::string pbcName;
        std::string algorithmName;
        // Make some symbolic names for the parameter combination.
        std::tie(numSettles, updateVelocities, calcVirial, pbcName, algorithmName) = GetParam();

        // Make a string that describes which parameter combination is
        // being tested, to help make failing tests comprehensible.
        std::string testDescription = formatString("while comparing CPU and GPU implementations with %d SETTLEs, %s, %svelocities and %scalculating the virial",
                                                   numSettles,
                                                   pbcName.c_str(),
                                                   updateVelocities ? "with " : "without ",
                                                   calcVirial ? "" : "not ");

        auto  testDataCpu = std::make_unique<SettleTestData>(numSettles);
        auto  testDataGpu = std::make_unique<SettleTestData>(numSettles);

        t_pbc pbc = pbcs_.at(pbcName);

        applySettle(testDataCpu.get(),
                    pbc,
                    updateVelocities,
                    calcVirial,
                    testDescription);

        applySettleGpu(testDataGpu.get(),
                       pbc,
                       updateVelocities,
                       calcVirial,
                       testDescription);

        FloatingPointTolerance tolerance = absoluteTolerance(0.0001);

        checkCoordinatesEqual(testDataGpu->xPrime_,
                              testDataCpu->xPrime_,
                              tolerance,
                              testDescription);

        if (updateVelocities)
        {
            checkCoordinatesEqual(testDataGpu->v_,
                                  testDataCpu->v_,
                                  tolerance,
                                  testDescription);
        }
        if (calcVirial)
        {
            checkTensorsEqual(testDataGpu->virial_,
                              testDataCpu->virial_,
                              tolerance,
                              testDescription);
        }
    }
}
#endif      // GMX_GPU == GMX_GPU_CUDA

// Scan the full Cartesian product of numbers of SETTLE interactions, whether or not update velocities
// or calculate the virial contribution, use or not to use PBC, use CPU or GPU-based version of SETTLE.
INSTANTIATE_TEST_CASE_P(WithParameters, SettleTest,
                            ::testing::Combine(::testing::Values(1, 2, 4, 5, 7, 10, 12, 15, 17),
                                                   ::testing::Bool(),
                                                   ::testing::Bool(),
                                                   ::testing::Values("PBCNone", "PBCXYZ"),
                                                   ::testing::ValuesIn(getAlgorithmsNames())
                                               ));
}  // namespace
}  // namespace test
}  // namespace gmx
