/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2020, by the GROMACS development team, led by
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
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include <string>

#include <gmock/gmock.h>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "pmetestcommon.h"
#include "testhardwarecontexts.h"

namespace gmx
{
namespace test
{
namespace
{

//! PME spline and spread code path being tested
enum class PmeSplineAndSpreadOptions
{
    SplineOnly,
    SpreadOnly,
    SplineAndSpreadUnified
};

/*! \brief Convenience typedef of input parameters - unit cell box, PME interpolation order, grid
 * dimensions, particle coordinates, particle charges
 * TODO: consider inclusion of local grid offsets/sizes or PME nodes counts to test the PME DD
 */
typedef std::tuple<Matrix3x3, int, IVec, CoordinatesVector, ChargesVector> SplineAndSpreadInputParameters;

/*! \brief Test fixture for testing both atom spline parameter computation and charge spreading.
 * These 2 stages of PME are tightly coupled in the code.
 */
class PmeSplineAndSpreadTest : public ::testing::TestWithParam<SplineAndSpreadInputParameters>
{
public:
    PmeSplineAndSpreadTest() = default;
    //! The test
    void runTest()
    {
        /* Getting the input */
        Matrix3x3         box;
        int               pmeOrder;
        IVec              gridSize;
        CoordinatesVector coordinates;
        ChargesVector     charges;

        std::tie(box, pmeOrder, gridSize, coordinates, charges) = GetParam();
        const size_t atomCount                                  = coordinates.size();

        /* Storing the input where it's needed */
        t_inputrec inputRec;
        inputRec.nkx         = gridSize[XX];
        inputRec.nky         = gridSize[YY];
        inputRec.nkz         = gridSize[ZZ];
        inputRec.pme_order   = pmeOrder;
        inputRec.coulombtype = eelPME;
        inputRec.epsilon_r   = 1.0;

        TestReferenceData refData;

        const std::map<PmeSplineAndSpreadOptions, std::string> optionsToTest = {
            { PmeSplineAndSpreadOptions::SplineAndSpreadUnified,
              "spline computation and charge spreading (fused)" },
            { PmeSplineAndSpreadOptions::SplineOnly, "spline computation" },
            { PmeSplineAndSpreadOptions::SpreadOnly, "charge spreading" }
        };

        // There is a subtle problem with multiple comparisons against same reference data:
        // The subsequent (GPU) spreading runs at one point didn't actually copy the output grid
        // into the proper buffer, but the reference data was already marked as checked
        // (hasBeenChecked_) by the CPU run, so nothing failed. For now we will manually track that
        // the count of the grid entries is the same on each run. This is just a hack for a single
        // specific output though. What would be much better TODO is to split different codepaths
        // into separate tests, while making them use the same reference files.
        bool   gridValuesSizeAssigned = false;
        size_t previousGridValuesSize;

        for (const auto& context : getPmeTestEnv()->getHardwareContexts())
        {
            CodePath   codePath = context->codePath();
            const bool supportedInput =
                    pmeSupportsInputForMode(*getPmeTestEnv()->hwinfo(), &inputRec, codePath);
            if (!supportedInput)
            {
                /* Testing the failure for the unsupported input */
                EXPECT_THROW_GMX(pmeInitWrapper(&inputRec, codePath, nullptr, nullptr, nullptr, box),
                                 NotImplementedError);
                continue;
            }

            for (const auto& option : optionsToTest)
            {
                /* Describing the test uniquely in case it fails */

                SCOPED_TRACE(formatString(
                        "Testing %s with %s %sfor PME grid size %d %d %d"
                        ", order %d, %zu atoms",
                        option.second.c_str(), codePathToString(codePath), context->description().c_str(),
                        gridSize[XX], gridSize[YY], gridSize[ZZ], pmeOrder, atomCount));

                /* Running the test */

                PmeSafePointer pmeSafe =
                        pmeInitWrapper(&inputRec, codePath, context->deviceContext(),
                                       context->deviceStream(), context->pmeGpuProgram(), box);
                std::unique_ptr<StatePropagatorDataGpu> stateGpu =
                        (codePath == CodePath::GPU)
                                ? makeStatePropagatorDataGpu(*pmeSafe.get(), context->deviceContext(),
                                                             context->deviceStream())
                                : nullptr;

                pmeInitAtoms(pmeSafe.get(), stateGpu.get(), codePath, coordinates, charges);

                const bool computeSplines =
                        (option.first == PmeSplineAndSpreadOptions::SplineOnly)
                        || (option.first == PmeSplineAndSpreadOptions::SplineAndSpreadUnified);
                const bool spreadCharges =
                        (option.first == PmeSplineAndSpreadOptions::SpreadOnly)
                        || (option.first == PmeSplineAndSpreadOptions::SplineAndSpreadUnified);

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
                        splineValuesChecker.checkSequence(splineValuesDim.begin(),
                                                          splineValuesDim.end(), dimString[i]);
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
                    rootChecker.checkSequence(gridLineIndices.begin(), gridLineIndices.end(),
                                              "Gridline indices");
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
};


/*! \brief Test for spline parameter computation and charge spreading. */
TEST_P(PmeSplineAndSpreadTest, ReproducesOutputs)
{
    EXPECT_NO_THROW_GMX(runTest());
}

/* Valid input instances */

//! A couple of valid inputs for boxes.
std::vector<Matrix3x3> const c_sampleBoxes{
    // normal box
    Matrix3x3{ { 8.0F, 0.0F, 0.0F, 0.0F, 3.4F, 0.0F, 0.0F, 0.0F, 2.0F } },
    // triclinic box
    Matrix3x3{ { 7.0F, 0.0F, 0.0F, 0.0F, 4.1F, 0.0F, 3.5F, 2.0F, 12.2F } },
};

//! A couple of valid inputs for grid sizes.
std::vector<IVec> const c_sampleGridSizes{ IVec{ 16, 12, 14 }, IVec{ 19, 17, 11 } };

//! Random charges
std::vector<real> const c_sampleChargesFull{ 4.95F, 3.11F, 3.97F, 1.08F, 2.09F, 1.1F,
                                             4.13F, 3.31F, 2.8F,  5.83F, 5.09F, 6.1F,
                                             2.86F, 0.24F, 5.76F, 5.19F, 0.72F };
//! 1 charge
auto const c_sampleCharges1 = ChargesVector(c_sampleChargesFull).subArray(0, 1);
//! 2 charges
auto const c_sampleCharges2 = ChargesVector(c_sampleChargesFull).subArray(1, 2);
//! 13 charges
auto const c_sampleCharges13 = ChargesVector(c_sampleChargesFull).subArray(3, 13);

//! Random coordinate vectors
CoordinatesVector const c_sampleCoordinatesFull{ { 5.59F, 1.37F, 0.95F },
                                                 {
                                                         16.0F, 1.02F, 0.22F // 2 box lengths in x
                                                 },
                                                 { 0.034F, 1.65F, 0.22F },
                                                 { 0.33F, 0.92F, 1.56F },
                                                 { 1.16F, 0.75F, 0.39F },
                                                 { 0.5F, 1.63F, 1.14F },
                                                 {
                                                         16.0001F, 1.52F, 1.19F // > 2 box lengths in x
                                                 },
                                                 {
                                                         1.43F, 1.1F, 4.1F // > 2 box lengths in z
                                                 },
                                                 {
                                                         -1.08F, 1.19F, 0.08F // negative x
                                                 },
                                                 { 1.6F, 0.93F, 0.53F },
                                                 {
                                                         1.32F, -1.48F, 0.16F // negative y
                                                 },
                                                 { 0.87F, 0.0F, 0.33F },
                                                 {
                                                         0.95F, 7.7F, -0.48F // > 2 box lengths in y, negative z
                                                 },
                                                 { 1.23F, 0.91F, 0.68F },
                                                 { 0.19F, 1.45F, 0.94F },
                                                 { 1.28F, 0.46F, 0.38F },
                                                 { 1.21F, 0.23F, 1.0F } };
//! 1 coordinate vector
CoordinatesVector const c_sampleCoordinates1(c_sampleCoordinatesFull.begin(),
                                             c_sampleCoordinatesFull.begin() + 1);
//! 2 coordinate vectors
CoordinatesVector const c_sampleCoordinates2(c_sampleCoordinatesFull.begin() + 1,
                                             c_sampleCoordinatesFull.begin() + 3);
//! 13 coordinate vectors
CoordinatesVector const c_sampleCoordinates13(c_sampleCoordinatesFull.begin() + 3,
                                              c_sampleCoordinatesFull.begin() + 16);

//! moved out from instantiantions for readability
auto c_inputBoxes = ::testing::ValuesIn(c_sampleBoxes);
//! moved out from instantiantions for readability
auto c_inputPmeOrders = ::testing::Range(3, 5 + 1);
//! moved out from instantiantions for readability
auto c_inputGridSizes = ::testing::ValuesIn(c_sampleGridSizes);

/*! \brief Instantiation of the test with valid input and 1 atom */
INSTANTIATE_TEST_CASE_P(SaneInput1,
                        PmeSplineAndSpreadTest,
                        ::testing::Combine(c_inputBoxes,
                                           c_inputPmeOrders,
                                           c_inputGridSizes,
                                           ::testing::Values(c_sampleCoordinates1),
                                           ::testing::Values(c_sampleCharges1)));

/*! \brief Instantiation of the test with valid input and 2 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput2,
                        PmeSplineAndSpreadTest,
                        ::testing::Combine(c_inputBoxes,
                                           c_inputPmeOrders,
                                           c_inputGridSizes,
                                           ::testing::Values(c_sampleCoordinates2),
                                           ::testing::Values(c_sampleCharges2)));
/*! \brief Instantiation of the test with valid input and 13 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput13,
                        PmeSplineAndSpreadTest,
                        ::testing::Combine(c_inputBoxes,
                                           c_inputPmeOrders,
                                           c_inputGridSizes,
                                           ::testing::Values(c_sampleCoordinates13),
                                           ::testing::Values(c_sampleCharges13)));
} // namespace
} // namespace test
} // namespace gmx
