/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
 * Implements PME spline computation and charge spreading tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include <string>

#include <gmock/gmock.h>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

#include "pmetestcommon.h"

namespace gmx
{
namespace test
{
namespace
{

//! A couple of valid inputs for grid sizes
std::vector<IVec> const c_inputGridSizes{ IVec{ 16, 12, 14 }, IVec{ 19, 17, 11 } };

//! PME spline and spread code path being tested
enum class SplineAndSpreadOptions
{
    SplineOnly,
    SpreadOnly,
    SplineAndSpreadUnified,
    Count
};

struct TestSystem
{
    CoordinatesVector coordinates;
    ChargesVector     charges;
};

const std::unordered_map<std::string, TestSystem> c_testSystems = {
    { "1 atom", { CoordinatesVector{ { 5.59F, 1.37F, 0.95F } }, ChargesVector{ 4.95F } } },
    { "2 atoms",
      { CoordinatesVector{
                // 2 box lengths in x
                { 16.0F, 1.02F, 0.22F },
                { 0.034F, 1.65F, 0.22F },
        },
        ChargesVector{ {
                3.11F,
                3.97F,
        } } } },
    { "13 atoms",
      { CoordinatesVector{
                { 0.33F, 0.92F, 1.56F },
                { 1.16F, 0.75F, 0.39F },
                { 0.5F, 1.63F, 1.14F },
                // > 2 box lengths in x
                { 16.0001F, 1.52F, 1.19F },
                // > 2 box lengths in z
                { 1.43F, 1.1F, 4.1F },
                // negative x
                { -1.08F, 1.19F, 0.08F },
                { 1.6F, 0.93F, 0.53F },
                // negative y
                { 1.32F, -1.48F, 0.16F },
                { 0.87F, 0.0F, 0.33F },
                // > 2 box lengths in y, negative z
                { 0.95F, 7.7F, -0.48F },
                { 1.23F, 0.91F, 0.68F },
                { 0.19F, 1.45F, 0.94F },
                { 1.28F, 0.46F, 0.38F },
        },
        ChargesVector{ 1.08F, 2.09F, 1.1F, 4.13F, 3.31F, 2.8F, 5.83F, 5.09F, 6.1F, 2.86F, 0.24F, 5.76F, 5.19F } } },
};

/* Valid input instances */

/*! \brief Convenience typedef of input parameters - unit cell box, PME interpolation order, grid
 * dimensions, particle coordinates, particle charges
 * TODO: consider inclusion of local grid offsets/sizes or PME nodes counts to test the PME DD
 */
typedef std::tuple<std::string, int, IVec, std::string> SplineAndSpreadInputParameters;

//! Help GoogleTest name our test cases
std::string nameOfTest(const testing::TestParamInfo<SplineAndSpreadInputParameters>& info)
{
    std::string testName = formatString(
            "box_%s_"
            "order_%d_"
            "grid_%d_%d_%d_"
            "system_%s",
            std::get<0>(info.param).c_str(),
            std::get<1>(info.param),
            std::get<2>(info.param)[XX],
            std::get<2>(info.param)[YY],
            std::get<2>(info.param)[ZZ],
            std::get<3>(info.param).c_str());

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");
    return testName;
}

const char* enumValueToString(SplineAndSpreadOptions enumValue)
{
    static constexpr gmx::EnumerationArray<SplineAndSpreadOptions, const char*> s_splineAndSpreadOptionsNames = {
        "spline computation",
        "charge spreading",
        "spline computation and charge spreading (fused)",
    };
    return s_splineAndSpreadOptionsNames[enumValue];
}

/*! \brief Test fixture for testing both atom spline parameter computation and charge spreading.
 * These 2 stages of PME are tightly coupled in the code.
 */
class SplineAndSpreadTest : public ::testing::TestWithParam<SplineAndSpreadInputParameters>
{
public:
    SplineAndSpreadTest() = default;

    static void SetUpTestSuite()
    {
        s_pmeTestHardwareContexts    = createPmeTestHardwareContextList();
        g_allowPmeWithSyclForTesting = true; // We support PmeSplineAndSpread with SYCL
    }

    static void TearDownTestSuite()
    {
        // Revert the value back.
        g_allowPmeWithSyclForTesting = false;
    }

    //! The test
    static void runTest()
    {
        /* Getting the input */
        int         pmeOrder;
        IVec        gridSize;
        std::string boxName, testSystemName;

        std::tie(boxName, pmeOrder, gridSize, testSystemName) = GetParam();
        Matrix3x3                box                          = c_inputBoxes.at(boxName);
        const CoordinatesVector& coordinates = c_testSystems.at(testSystemName).coordinates;
        const ChargesVector&     charges     = c_testSystems.at(testSystemName).charges;
        const size_t             atomCount   = coordinates.size();

        /* Storing the input where it's needed */
        t_inputrec inputRec;
        inputRec.nkx         = gridSize[XX];
        inputRec.nky         = gridSize[YY];
        inputRec.nkz         = gridSize[ZZ];
        inputRec.pme_order   = pmeOrder;
        inputRec.coulombtype = CoulombInteractionType::Pme;
        inputRec.epsilon_r   = 1.0;

        // There is a subtle problem with multiple comparisons against same reference data:
        // The subsequent (GPU) spreading runs at one point didn't actually copy the output grid
        // into the proper buffer, but the reference data was already marked as checked
        // (hasBeenChecked_) by the CPU run, so nothing failed. For now we will manually track that
        // the count of the grid entries is the same on each run. This is just a hack for a single
        // specific output though. What would be much better TODO is to split different codepaths
        // into separate tests, while making them use the same reference files.
        bool   gridValuesSizeAssigned = false;
        size_t previousGridValuesSize;

        TestReferenceData refData;
        for (const auto& pmeTestHardwareContext : s_pmeTestHardwareContexts)
        {
            pmeTestHardwareContext->activate();
            CodePath   codePath       = pmeTestHardwareContext->codePath();
            const bool supportedInput = pmeSupportsInputForMode(
                    *getTestHardwareEnvironment()->hwinfo(), &inputRec, codePath);
            if (!supportedInput)
            {
                /* Testing the failure for the unsupported input */
                EXPECT_THROW_GMX(pmeInitWrapper(&inputRec, codePath, nullptr, nullptr, nullptr, box),
                                 NotImplementedError);
                continue;
            }

            for (SplineAndSpreadOptions option : EnumerationWrapper<SplineAndSpreadOptions>())
            {
                /* Describing the test uniquely in case it fails */

                SCOPED_TRACE(
                        formatString("Testing %s on %s for PME grid size %d %d %d"
                                     ", order %d, %zu atoms",
                                     enumValueToString(option),
                                     pmeTestHardwareContext->description().c_str(),
                                     gridSize[XX],
                                     gridSize[YY],
                                     gridSize[ZZ],
                                     pmeOrder,
                                     atomCount));

                /* Running the test */

                PmeSafePointer                          pmeSafe = pmeInitWrapper(&inputRec,
                                                        codePath,
                                                        pmeTestHardwareContext->deviceContext(),
                                                        pmeTestHardwareContext->deviceStream(),
                                                        pmeTestHardwareContext->pmeGpuProgram(),
                                                        box);
                std::unique_ptr<StatePropagatorDataGpu> stateGpu =
                        (codePath == CodePath::GPU)
                                ? makeStatePropagatorDataGpu(*pmeSafe.get(),
                                                             pmeTestHardwareContext->deviceContext(),
                                                             pmeTestHardwareContext->deviceStream())
                                : nullptr;

                pmeInitAtoms(pmeSafe.get(), stateGpu.get(), codePath, coordinates, charges);

                const bool computeSplines = (option == SplineAndSpreadOptions::SplineOnly)
                                            || (option == SplineAndSpreadOptions::SplineAndSpreadUnified);
                const bool spreadCharges = (option == SplineAndSpreadOptions::SpreadOnly)
                                           || (option == SplineAndSpreadOptions::SplineAndSpreadUnified);

                if (!computeSplines)
                {
                    // Here we should set up the results of the spline computation so that the spread can run.
                    // What is lazy and works is running the separate spline so that it will set it up for us:
                    pmePerformSplineAndSpread(pmeSafe.get(), codePath, true, false);
                    // We know that it is tested in another iteration.
                    // TODO: Clean alternative: read and set the reference gridline indices, spline params
                }

                pmePerformSplineAndSpread(pmeSafe.get(), codePath, computeSplines, spreadCharges);
                pmeFinalizeTest(pmeSafe.get(), codePath);

                /* Outputs correctness check */
                /* All tolerances were picked empirically for single precision on CPU */

                TestReferenceChecker rootChecker(refData.rootChecker());

                const auto maxGridSize = std::max(std::max(gridSize[XX], gridSize[YY]), gridSize[ZZ]);
                const auto ulpToleranceSplineValues = 4 * (pmeOrder - 2) * maxGridSize;
                /* 4 is a modest estimate for amount of operations; (pmeOrder - 2) is a number of iterations;
                 * maxGridSize is inverse of the smallest positive fractional coordinate (which are interpolated by the splines).
                 */

                if (computeSplines)
                {
                    const char* dimString[] = { "X", "Y", "Z" };

                    /* Spline values */
                    SCOPED_TRACE(formatString("Testing spline values with tolerance of %d",
                                              ulpToleranceSplineValues));
                    TestReferenceChecker splineValuesChecker(
                            rootChecker.checkCompound("Splines", "Values"));
                    splineValuesChecker.setDefaultTolerance(
                            relativeToleranceAsUlp(1.0, ulpToleranceSplineValues));
                    for (int i = 0; i < DIM; i++)
                    {
                        auto splineValuesDim =
                                pmeGetSplineData(pmeSafe.get(), codePath, PmeSplineDataType::Values, i);
                        splineValuesChecker.checkSequence(
                                splineValuesDim.begin(), splineValuesDim.end(), dimString[i]);
                    }

                    /* Spline derivatives */
                    const auto ulpToleranceSplineDerivatives = 4 * ulpToleranceSplineValues;
                    /* 4 is just a wild guess since the derivatives are deltas of neighbor spline values which could differ greatly */
                    SCOPED_TRACE(formatString("Testing spline derivatives with tolerance of %d",
                                              ulpToleranceSplineDerivatives));
                    TestReferenceChecker splineDerivativesChecker(
                            rootChecker.checkCompound("Splines", "Derivatives"));
                    splineDerivativesChecker.setDefaultTolerance(
                            relativeToleranceAsUlp(1.0, ulpToleranceSplineDerivatives));
                    for (int i = 0; i < DIM; i++)
                    {
                        auto splineDerivativesDim = pmeGetSplineData(
                                pmeSafe.get(), codePath, PmeSplineDataType::Derivatives, i);
                        splineDerivativesChecker.checkSequence(
                                splineDerivativesDim.begin(), splineDerivativesDim.end(), dimString[i]);
                    }

                    /* Particle gridline indices */
                    auto gridLineIndices = pmeGetGridlineIndices(pmeSafe.get(), codePath);
                    rootChecker.checkSequence(
                            gridLineIndices.begin(), gridLineIndices.end(), "Gridline indices");
                }

                if (spreadCharges)
                {
                    /* The wrapped grid */
                    SparseRealGridValuesOutput nonZeroGridValues = pmeGetRealGrid(pmeSafe.get(), codePath);
                    TestReferenceChecker gridValuesChecker(
                            rootChecker.checkCompound("NonZeroGridValues", "RealSpaceGrid"));
                    const auto ulpToleranceGrid =
                            2 * ulpToleranceSplineValues
                            * static_cast<int>(ceil(sqrt(static_cast<real>(atomCount))));
                    /* 2 is empiric; sqrt(atomCount) assumes all the input charges may spread onto the same cell */
                    SCOPED_TRACE(formatString("Testing grid values with tolerance of %d", ulpToleranceGrid));
                    if (!gridValuesSizeAssigned)
                    {
                        previousGridValuesSize = nonZeroGridValues.size();
                        gridValuesSizeAssigned = true;
                    }
                    else
                    {
                        EXPECT_EQ(previousGridValuesSize, nonZeroGridValues.size());
                    }

                    gridValuesChecker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpToleranceGrid));
                    for (const auto& point : nonZeroGridValues)
                    {
                        gridValuesChecker.checkReal(point.second, point.first.c_str());
                    }
                }
            }
        }
    }

    static std::vector<std::unique_ptr<PmeTestHardwareContext>> s_pmeTestHardwareContexts;
};

std::vector<std::unique_ptr<PmeTestHardwareContext>> SplineAndSpreadTest::s_pmeTestHardwareContexts;


/*! \brief Test for spline parameter computation and charge spreading. */
TEST_P(SplineAndSpreadTest, WorksWith)
{
    EXPECT_NO_THROW_GMX(runTest());
}

//! Moved out from instantiations for readability
const auto c_inputBoxNames = ::testing::Values("rect", "tric");
//! Moved out from instantiations for readability
const auto c_inputGridNames = ::testing::Values("first", "second");
//! Moved out from instantiations for readability
const auto c_inputTestSystemNames = ::testing::Values("1 atom", "2 atoms", "13 atoms");

INSTANTIATE_TEST_SUITE_P(Pme,
                         SplineAndSpreadTest,
                         ::testing::Combine(c_inputBoxNames,
                                            ::testing::ValuesIn(c_inputPmeOrders),
                                            ::testing::ValuesIn(c_inputGridSizes),
                                            c_inputTestSystemNames),
                         nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
