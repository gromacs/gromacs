/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdlib/tests/watersystem.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/hardware_test_fixture.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "settletestdata.h"

#if GMX_GPU && !GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/device_context.h"
#    include "gromacs/gpu_utils/device_stream.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/mdlib/settle_gpu.h"
#endif

namespace gmx
{
namespace test
{
namespace
{

void applySettleCpu(SettleTestData*    testData,
                    const t_pbc&       pbc,
                    const bool         updateVelocities,
                    const bool         calcVirial,
                    const std::string& testDescription)
{
    SettleData settled(testData->mtop_);

    settled.setConstraints(testData->idef_->il[InteractionFunction::SETTLE],
                           testData->numAtoms_,
                           testData->masses_,
                           testData->inverseMasses_);

    bool errorOccured;
    int  numThreads  = 1;
    int  threadIndex = 0;
    csettle(settled,
            numThreads,
            threadIndex,
            &pbc,
            testData->x_.arrayRefWithPadding(),
            testData->xPrime_.arrayRefWithPadding(),
            testData->reciprocalTimeStep_,
            updateVelocities ? testData->v_.arrayRefWithPadding() : ArrayRefWithPadding<RVec>(),
            calcVirial,
            testData->virial_,
            &errorOccured);
    EXPECT_FALSE(errorOccured) << testDescription;
}

#if GMX_GPU && !GMX_GPU_OPENCL

void applySettleGpu(const DeviceContext& deviceContext,
                    const DeviceStream&  deviceStream,
                    SettleTestData*      testData,
                    const t_pbc&         pbc,
                    const bool           updateVelocities,
                    const bool           calcVirial,
                    const std::string& /* testDescription */)
{
    deviceContext.activate();

    auto settleGpu = std::make_unique<SettleGpu>(testData->mtop_, deviceContext, deviceStream);

    settleGpu->set(*testData->idef_);
    PbcAiuc pbcAiuc;
    setPbcAiuc(pbc.ndim_ePBC, pbc.box, &pbcAiuc);

    int numAtoms = testData->numAtoms_;

    DeviceBuffer<Float3> d_x, d_xp, d_v;

    Float3* h_x  = gmx::asGenericFloat3Pointer(testData->x_);
    Float3* h_xp = gmx::asGenericFloat3Pointer(testData->xPrime_);
    Float3* h_v  = gmx::asGenericFloat3Pointer(testData->v_);

    allocateDeviceBuffer(&d_x, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_xp, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_v, numAtoms, deviceContext);

    copyToDeviceBuffer(&d_x, h_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, h_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyToDeviceBuffer(&d_v, h_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    }
    settleGpu->apply(
            d_x, d_xp, updateVelocities, d_v, testData->reciprocalTimeStep_, calcVirial, testData->virial_, pbcAiuc);

    copyFromDeviceBuffer(h_xp, &d_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyFromDeviceBuffer(h_v, &d_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    }

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
}

#endif // GMX_GPU && !GMX_GPU_OPENCL

//! PBC types for testing
enum class PbcTestType : int
{
    None,
    XYZ,
    Count
};

static constexpr real dOHTest = 0.10;
static constexpr real dHHTest = 0.15;
static constexpr real mOTest  = 10;
static constexpr real mHTest  = 2;

//! Static const array of test PBCs, initialized once for all tests
static const gmx::EnumerationArray<PbcTestType, t_pbc> sc_testPbcs = []()
{
    gmx::EnumerationArray<PbcTestType, t_pbc> pbcs;
    t_pbc                                     pbc;

    // Infinitely small box (no PBC)
    matrix boxNone = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    set_pbc(&pbc, PbcType::No, boxNone);
    pbcs[PbcTestType::None] = pbc;

    // Rectangular box (XYZ PBC)
    matrix boxXyz = { { real(1.86206), 0, 0 }, { 0, real(1.86206), 0 }, { 0, 0, real(1.86206) } };
    set_pbc(&pbc, PbcType::Xyz, boxXyz);
    pbcs[PbcTestType::XYZ] = pbc;

    return pbcs;
}();

//! Names for PBC test types
static const gmx::EnumerationArray<PbcTestType, const char*> sc_pbcTestTypeNames = { { "PBCNone",
                                                                                       "PBCXYZ" } };

std::unique_ptr<gmx_mtop_t> mtopTwoIdenticalMoltypes()
{
    t_iparams ip1;
    ip1.settle.doh = dOHTest;
    ip1.settle.dhh = dHHTest;

    std::unique_ptr<gmx_mtop_t> mtop = std::make_unique<gmx_mtop_t>();
    mtop->ffparams.iparams.push_back(ip1);
    mtop->ffparams.iparams.push_back(ip1);

    const int nral = NRAL(InteractionFunction::SETTLE);

    mtop->moltype.resize(2);
    for (int m = 0; m < 2; m++)
    {
        gmx_moltype_t& molt = mtop->moltype[m];

        molt.atoms.nr = nral;
        snew(molt.atoms.atom, molt.atoms.nr);
        molt.atoms.atom[0].m = mOTest;
        molt.atoms.atom[1].m = mHTest;
        molt.atoms.atom[2].m = mHTest;

        molt.ilist[InteractionFunction::SETTLE].iatoms.push_back(0);
        for (int a = 0; a < NRAL(InteractionFunction::SETTLE); a++)
        {
            molt.ilist[InteractionFunction::SETTLE].iatoms.push_back(a);
        }
    }

    return mtop;
}

TEST(Settle, MultipleIdenticalSettlesWork)
{
    std::unique_ptr<gmx_mtop_t> mtop = mtopTwoIdenticalMoltypes();

    SettleWaterTopology settleTop = getSettleTopologyData(*mtop);

    EXPECT_EQ(settleTop.dOH, dOHTest);
    EXPECT_EQ(settleTop.dHH, dHHTest);
    EXPECT_EQ(settleTop.mO, mOTest);
    EXPECT_EQ(settleTop.mH, mHTest);
}

TEST(Settle, MultipleDifferentSettlesThrow)
{
    std::unique_ptr<gmx_mtop_t> mtop = mtopTwoIdenticalMoltypes();
    mtop->moltype[1].atoms.atom[0].m *= 0.5;

    EXPECT_THROW(getSettleTopologyData(*mtop), InvalidInputError);
}

//! Input configuration for SETTLE tests (defines the physical scenario being tested)
using SettleInputConfig = std::tuple<int, bool, bool, PbcTestType>;

/*! Hardware test helper for SETTLE
 *
 * \todo There are no execution modes - tests don't vary SETTLE
 * implementation, but they should test SIMD vs no SIMD here. */
using SettleTestHelper = HardwareAndExecutionTestHelper<SettleInputConfig, std::tuple<>>;

//! Formatters for parameters in the config info
static const auto sc_configInfoFormatters =
        std::make_tuple([](int n) { return formatString("%dsettles", n); },
                        [](bool b) { return b ? "velocities" : "novelocities"; },
                        [](bool b) { return b ? "virial" : "novirial"; },
                        [](PbcTestType pbc) { return sc_pbcTestTypeNames[pbc]; });
//! Formatters for parameters in the execution mode (currently empty)
static const auto sc_executionModeFormatters = std::make_tuple();

//! Helper object to name tests using all parameters
static const NameOfTestFromTuple<SettleTestHelper::DynamicParameters> sc_testNamer =
        SettleTestHelper::testNamer(sc_configInfoFormatters, sc_executionModeFormatters);

/*! \brief Sets of parameters on which to run the tests.
 *
 * Each entry is a tuple of (numSettles, updateVelocities, calcVirial, pbcType).
 */
const std::array<SettleInputConfig, 13> sc_settleConfigs = { {
        { 1, false, false, PbcTestType::XYZ },   // 1 water molecule
        { 2, false, false, PbcTestType::XYZ },   // 2 water molecules
        { 4, false, false, PbcTestType::XYZ },   // 4 water molecules
        { 5, false, false, PbcTestType::XYZ },   // 5 water molecules
        { 6, false, false, PbcTestType::XYZ },   // 6 water molecules
        { 10, false, false, PbcTestType::XYZ },  // 10 water molecules
        { 12, false, false, PbcTestType::XYZ },  // 12 water molecules
        { 15, false, false, PbcTestType::XYZ },  // 15 water molecules
        { 17, true, false, PbcTestType::XYZ },   // Update velocities
        { 17, false, true, PbcTestType::XYZ },   // Compute virial
        { 17, false, false, PbcTestType::None }, // No periodic boundary
        { 17, true, true, PbcTestType::None },   // Update velocities, compute virial, without PBC
        { 17, true, true, PbcTestType::XYZ }     // Update velocities, compute virial, with PBC
} };

/*! \brief Test fixture for testing SETTLE.
 */
class SettleTest : public HardwareTestFixture<SettleTestHelper>
{
public:
    SettleTest() : HardwareTestFixture(sc_configInfoFormatters) {}

    /*! \brief Check if the final interatomic distances are equal to target set by constraints.
     *
     * \param[in]  numSettles        Number of water molecules in the tested system.
     * \param[in]  tolerance         Tolerance to compare floating point numbers.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    static void checkConstrainsSatisfied(const int                    numSettles,
                                         const FloatingPointTolerance tolerance,
                                         const SettleTestData&        testData)
    {
        for (int i = 0; i < numSettles; ++i)
        {
            const gmx::RVec& positionO  = testData.xPrime_[i * testData.atomsPerSettle_ + 0];
            const gmx::RVec& positionH1 = testData.xPrime_[i * testData.atomsPerSettle_ + 1];
            const gmx::RVec& positionH2 = testData.xPrime_[i * testData.atomsPerSettle_ + 2];

            real dOH = testData.dOH_;
            real dHH = testData.dHH_;

            EXPECT_REAL_EQ_TOL(dOH * dOH, distance2(positionO, positionH1), tolerance)
                    << formatString("for water %d. ", i);
            EXPECT_REAL_EQ_TOL(dOH * dOH, distance2(positionO, positionH2), tolerance)
                    << formatString("for water %d. ", i);
            EXPECT_REAL_EQ_TOL(dHH * dHH, distance2(positionH1, positionH2), tolerance)
                    << formatString("for water %d. ", i);
        }
    }

    /*! \brief Check if the virial was updated and symmetric.
     *
     * The two tests on virial are:
     * 1. If it was updated in case calcVirial is true.
     * 2. If it is symmetrical.
     *
     * \param[in]  calcVirial        If the virial is computed.
     * \param[in]  tolerance         Tolerance to compare floating point numbers.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    static void checkVirialSymmetric(const bool                   calcVirial,
                                     const FloatingPointTolerance tolerance,
                                     const SettleTestData&        testData)
    {
        for (int d = 0; d < DIM; ++d)
        {
            for (int dd = 0; dd < DIM; ++dd)
            {

                EXPECT_TRUE(calcVirial == (0. != testData.virial_[d][dd]))
                        << formatString("for virial component[%d][%d]. ", d, dd);

                if (calcVirial)
                {
                    EXPECT_REAL_EQ_TOL(testData.virial_[d][dd], testData.virial_[dd][d], tolerance)
                            << formatString("Virial is not symmetrical for [%d][%d]. ", d, dd);
                }
            }
        }
    }

    /*! \brief Check if the final positions correspond to reference values.
     *
     * \param[in]  numSettles        Number of water molecules in the tested system.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    void checkFinalPositions(const int numSettles, const SettleTestData& testData)
    {
        TestReferenceChecker finalCoordinatesRef(
                checker().checkSequenceCompound("FinalCoordinates", numSettles));
        for (int i = 0; i < numSettles; ++i)
        {
            TestReferenceChecker settlerRef(finalCoordinatesRef.checkCompound("Settler", nullptr));
            TestReferenceChecker atomsRef(
                    settlerRef.checkSequenceCompound("Atoms", testData.atomsPerSettle_));
            for (int j = 0; j < testData.atomsPerSettle_; ++j)
            {
                const gmx::RVec&     xPrime = testData.xPrime_[testData.atomsPerSettle_ * i + j];
                TestReferenceChecker xPrimeRef(atomsRef.checkCompound("Atom", nullptr));
                xPrimeRef.checkReal(xPrime[XX], "XX");
                xPrimeRef.checkReal(xPrime[YY], "YY");
                xPrimeRef.checkReal(xPrime[ZZ], "ZZ");
            }
        }
    }

    /*! \brief Check if the final velocities correspond to reference values.
     *
     * \param[in]  numSettles        Number of water molecules in the tested system.
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    void checkFinalVelocities(const int numSettles, const SettleTestData& testData)
    {
        TestReferenceChecker finalCoordinatesRef(
                checker().checkSequenceCompound("FinalVelocities", numSettles));
        for (int i = 0; i < numSettles; ++i)
        {
            TestReferenceChecker settlerRef(finalCoordinatesRef.checkCompound("Settler", nullptr));
            TestReferenceChecker atomsRef(
                    settlerRef.checkSequenceCompound("Atoms", testData.atomsPerSettle_));
            for (int j = 0; j < testData.atomsPerSettle_; ++j)
            {
                const gmx::RVec&     v = testData.v_[testData.atomsPerSettle_ * i + j];
                TestReferenceChecker vRef(atomsRef.checkCompound("Atom", nullptr));
                vRef.checkReal(v[XX], "XX");
                vRef.checkReal(v[YY], "YY");
                vRef.checkReal(v[ZZ], "ZZ");
            }
        }
    }

    /*! \brief Check if the computed virial correspond to reference values.
     *
     * \param[in]  testData          An object, containing all the data structures needed by SETTLE.
     */
    void checkVirial(const SettleTestData& testData)
    {
        const tensor&        virial = testData.virial_;
        TestReferenceChecker virialRef(checker().checkCompound("Virial", nullptr));

        // TODO: Is it worth it to make this in a loop??
        virialRef.checkReal(virial[XX][XX], "XX");
        virialRef.checkReal(virial[XX][YY], "XY");
        virialRef.checkReal(virial[XX][ZZ], "XZ");
        virialRef.checkReal(virial[YY][XX], "YX");
        virialRef.checkReal(virial[YY][YY], "YY");
        virialRef.checkReal(virial[YY][ZZ], "YZ");
        virialRef.checkReal(virial[ZZ][XX], "ZX");
        virialRef.checkReal(virial[ZZ][YY], "ZY");
        virialRef.checkReal(virial[ZZ][ZZ], "ZZ");
    }
};

TEST_P(SettleTest, SatisfiesConstraints)
{
    // Extract parameters (skip hardware context with _)
    auto [numSettles, updateVelocities, calcVirial, pbcType, _] = GetParam();

    // Extra indentation for reviewer convenience
    {
        const t_pbc& pbc = sc_testPbcs[pbcType];

        // Make a string that describes which parameter combination is being tested
        std::string testDescription = formatString(
                "Testing %s with %d SETTLEs, %s, %svelocities and %scalculating the virial.",
                hardwareContext()->description().c_str(),
                numSettles,
                sc_pbcTestTypeNames[pbcType],
                updateVelocities ? "with " : "without ",
                calcVirial ? "" : "not ");

        SCOPED_TRACE(testDescription);

        auto testData = std::make_unique<SettleTestData>(numSettles);

        ASSERT_LE(numSettles, testData->xPrime_.size() / testData->atomsPerSettle_)
                << "cannot test that many SETTLEs. " << testDescription;

        // Apply SETTLE based on hardware context
        if (isGpuTest())
        {
#if GMX_GPU && !GMX_GPU_OPENCL
            activateHardware();
            applySettleGpu(
                    *deviceContext(), *deviceStream(), testData.get(), pbc, updateVelocities, calcVirial, testDescription);
#else
            GMX_THROW(gmx::InternalError("GPU hardware context with no test code"));
#endif
        }
        else
        {
            applySettleCpu(testData.get(), pbc, updateVelocities, calcVirial, testDescription);
        }

        // The necessary tolerances for the test to pass were determined
        // empirically. This isn't nice, but the required behavior that
        // SETTLE produces constrained coordinates consistent with
        // sensible sampling needs to be tested at a much higher level.
        // TODO: Re-evaluate the tolerances.
        real dOH = testData->dOH_;
        FloatingPointTolerance tolerance = relativeToleranceAsPrecisionDependentUlp(dOH * dOH, 80, 380);
        FloatingPointTolerance toleranceVirial = absoluteTolerance(0.000001);

        FloatingPointTolerance tolerancePositions  = absoluteTolerance(0.000001);
        FloatingPointTolerance toleranceVelocities = absoluteTolerance(0.0001);

        checkConstrainsSatisfied(numSettles, tolerance, *testData);
        checkVirialSymmetric(calcVirial, toleranceVirial, *testData);

        checker().setDefaultTolerance(tolerancePositions);
        checkFinalPositions(numSettles, *testData);

        if (updateVelocities)
        {
            checker().setDefaultTolerance(toleranceVelocities);
            checkFinalVelocities(numSettles, *testData);
        }

        if (calcVirial)
        {
            checker().setDefaultTolerance(toleranceVirial);
            checkVirial(*testData);
        }
    }
}

// Run test on pre-determined set of combinations for test parameters, which include the numbers of SETTLEs (water
// molecules), whether or not velocities are updated and virial contribution is computed, was the PBC enabled.
// Each test runs as a separate GoogleTest case for each hardware configuration.
INSTANTIATE_TEST_SUITE_P(AllHardware,
                         SettleTest,
                         ::testing::ConvertGenerator(
                                 ::testing::Combine(::testing::ValuesIn(sc_settleConfigs),
                                                    ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                            GpuConfigurationCapabilities::Update))),
                                 flattenTupleWithHardwareContext<SettleInputConfig>()),
                         sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
