/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

/*! \brief Convenience typedef of input parameters - unit cell box, PME interpolation order, grid dimensions,
 * particle coordinates, particle charges
 * TODO: consider inclusion of local grid offsets/sizes or PME nodes counts to test the PME DD
 */
typedef std::tuple<Matrix3x3, int, IVec, CoordinatesVector, ChargesVector> SplineAndSpreadInputParameters;

/*! \brief Test fixture for testing both atom spline parameter computation and charge spreading.
 * These 2 stages of PME are tightly coupled in the code.
 */
class PmeSplineAndSpreadTest : public ::testing::TestWithParam<SplineAndSpreadInputParameters>
{
    public:
        //! Default constructor
        PmeSplineAndSpreadTest() = default;
        //! The test
        void runTest()
        {
            /* Getting the input */
            Matrix3x3               box;
            int                     pmeOrder;
            IVec                    gridSize;
            CoordinatesVector       coordinates;
            ChargesVector           charges;

            std::tie(box, pmeOrder, gridSize, coordinates, charges) = GetParam();
            const size_t atomCount = coordinates.size();

            /* Storing the input where it's needed */
            t_inputrec inputRec;
            inputRec.nkx         = gridSize[XX];
            inputRec.nky         = gridSize[YY];
            inputRec.nkz         = gridSize[ZZ];
            inputRec.pme_order   = pmeOrder;
            inputRec.coulombtype = eelPME;
            inputRec.epsilon_r   = 1.0;

            TestReferenceData                                      refData;

            const std::map<CodePath, std::string>                  modesToTest   = {{CodePath::CPU, "CPU"},
                                                                                    {CodePath::CUDA, "CUDA"}};

            const std::map<PmeSplineAndSpreadOptions, std::string> optionsToTest = {{PmeSplineAndSpreadOptions::SplineAndSpreadUnified, "spline computation and charge spreading (fused)"},
                                                                                    {PmeSplineAndSpreadOptions::SplineOnly, "spline computation"},
                                                                                    {PmeSplineAndSpreadOptions::SpreadOnly, "charge spreading"}};

            // There is a subtle problem with multiple comparisons against same reference data:
            // The subsequent (GPU) spreading runs at one point didn't actually copy the output grid into the proper buffer,
            // but the reference data was already marked as checked (hasBeenChecked_) by the CPU run, so nothing failed.
            // For now we will manually track that the count of the grid entries is the same on each run.
            // This is just a hack for a single specific output though.
            // What would be much better TODO is to split different codepaths into separate tests,
            // while making them use the same reference files.
            bool   gridValuesSizeAssigned = false;
            size_t previousGridValuesSize;

            for (const auto &mode : modesToTest)
            {
                const bool supportedInput = pmeSupportsInputForMode(&inputRec, mode.first);
                if (!supportedInput)
                {
                    /* Testing the failure for the unsupported input */
                    EXPECT_THROW(pmeInitAtoms(&inputRec, mode.first, nullptr, coordinates, charges, box), NotImplementedError);
                    continue;
                }

                const auto contextsToTest = pmeEnv->getHardwareContexts(mode.first);
                for (const auto &context : contextsToTest)
                {
                    for (const auto &option : optionsToTest)
                    {
                        /* Describing the test uniquely in case it fails */

                        SCOPED_TRACE(formatString("Testing %s with %s %sfor PME grid size %d %d %d"
                                                  ", order %d, %zu atoms",
                                                  option.second.c_str(), mode.second.c_str(),
                                                  context.getDescription().c_str(),
                                                  gridSize[XX], gridSize[YY], gridSize[ZZ],
                                                  pmeOrder,
                                                  atomCount));

                        /* Running the test */

                        PmeSafePointer pmeSafe = pmeInitAtoms(&inputRec, mode.first, context.getDeviceInfo(), coordinates, charges, box);

                        const bool     computeSplines = (option.first == PmeSplineAndSpreadOptions::SplineOnly) || (option.first == PmeSplineAndSpreadOptions::SplineAndSpreadUnified);
                        const bool     spreadCharges  = (option.first == PmeSplineAndSpreadOptions::SpreadOnly) || (option.first == PmeSplineAndSpreadOptions::SplineAndSpreadUnified);

                        if (!computeSplines)
                        {
                            // Here we should set up the results of the spline computation so that the spread can run.
                            // What is lazy and works is running the separate spline so that it will set it up for us:
                            pmePerformSplineAndSpread(pmeSafe.get(), mode.first, true, false);
                            // We know that it is tested in another iteration.
                            // TODO: Clean alternative: read and set the reference gridline indices, spline params
                        }

                        pmePerformSplineAndSpread(pmeSafe.get(), mode.first, computeSplines, spreadCharges);
                        pmeFinalizeTest(pmeSafe.get(), mode.first);

                        /* Outputs correctness check */
                        /* All tolerances were picked empirically for single precision on CPU */

                        TestReferenceChecker rootChecker(refData.rootChecker());

                        const auto           maxGridSize              = std::max(std::max(gridSize[XX], gridSize[YY]), gridSize[ZZ]);
                        const auto           ulpToleranceSplineValues = 2 * (pmeOrder - 2) * maxGridSize;
                        /* 2 is empiric, the rest follows from the amount of operations */

                        if (computeSplines)
                        {
                            const char *dimString[] = { "X", "Y", "Z" };

                            /* Spline values */
                            SCOPED_TRACE(formatString("Testing spline values with tolerance of %ld", ulpToleranceSplineValues));
                            TestReferenceChecker splineValuesChecker(rootChecker.checkCompound("Splines", "Values"));
                            splineValuesChecker.setDefaultTolerance(getSplineTolerance(ulpToleranceSplineValues));
                            for (int i = 0; i < DIM; i++)
                            {
                                auto splineValuesDim = pmeGetSplineData(pmeSafe.get(), mode.first, PmeSplineDataType::Values, i);
                                splineValuesChecker.checkSequence(splineValuesDim.begin(), splineValuesDim.end(), dimString[i]);
                            }

                            /* Spline derivatives */
                            const auto ulpToleranceSplineDerivatives = 4 * ulpToleranceSplineValues;
                            /* 4 is just a wild guess since the derivatives are deltas of neighbor spline values which could differ greatly */
                            SCOPED_TRACE(formatString("Testing spline derivatives with tolerance of %ld", ulpToleranceSplineDerivatives));
                            TestReferenceChecker splineDerivativesChecker(rootChecker.checkCompound("Splines", "Derivatives"));
                            splineDerivativesChecker.setDefaultTolerance(getSplineTolerance(ulpToleranceSplineDerivatives));
                            for (int i = 0; i < DIM; i++)
                            {
                                auto splineDerivativesDim = pmeGetSplineData(pmeSafe.get(), mode.first, PmeSplineDataType::Derivatives, i);
                                splineDerivativesChecker.checkSequence(splineDerivativesDim.begin(), splineDerivativesDim.end(), dimString[i]);
                            }

                            /* Particle gridline indices */
                            auto gridLineIndices = pmeGetGridlineIndices(pmeSafe.get(), mode.first);
                            rootChecker.checkSequence(gridLineIndices.begin(), gridLineIndices.end(), "Gridline indices");
                        }

                        if (spreadCharges)
                        {
                            /* The wrapped grid */
                            SparseRealGridValuesOutput nonZeroGridValues = pmeGetRealGrid(pmeSafe.get(), mode.first);
                            TestReferenceChecker       gridValuesChecker(rootChecker.checkCompound("NonZeroGridValues", "RealSpaceGrid"));
                            const auto                 ulpToleranceGrid = 2 * ulpToleranceSplineValues * (int)(ceil(sqrt(atomCount)));
                            /* 2 is empiric; sqrt(atomCount) assumes all the input charges may spread onto the same cell */
                            SCOPED_TRACE(formatString("Testing grid values with tolerance of %ld", ulpToleranceGrid));
                            if (!gridValuesSizeAssigned)
                            {
                                previousGridValuesSize = nonZeroGridValues.size();
                                gridValuesSizeAssigned = true;
                            }
                            else
                            {
                                EXPECT_EQ(previousGridValuesSize, nonZeroGridValues.size());
                            }

                            gridValuesChecker.setDefaultTolerance(getSplineTolerance(ulpToleranceGrid));
                            for (const auto &point : nonZeroGridValues)
                            {
                                gridValuesChecker.checkReal(point.second, point.first.c_str());
                            }
                        }
                    }
                }
            }
        }
};


/*! \brief Test for spline parameter computation and charge spreading. */
TEST_P(PmeSplineAndSpreadTest, ReproducesOutputs)
{
    EXPECT_NO_THROW(runTest());
}

/* Valid input instances */

//! A couple of valid inputs for boxes.
static std::vector<Matrix3x3> const c_sampleBoxes
{
    // normal box
    Matrix3x3 {{
                   8.0f, 0.0f, 0.0f,
                   0.0f, 3.4f, 0.0f,
                   0.0f, 0.0f, 2.0f
               }},
    // triclinic box
    Matrix3x3 {{
                   7.0f, 0.0f, 0.0f,
                   0.0f, 4.1f, 0.0f,
                   3.5f, 2.0f, 12.2f
               }},
};

//! A couple of valid inputs for grid sizes.
static std::vector<IVec> const c_sampleGridSizes
{
    IVec {
        16, 12, 14
    },
    IVec {
        19, 17, 11
    }
};

//! Random charges
static std::vector<real> const c_sampleChargesFull
{
    4.95f, 3.11f, 3.97f, 1.08f, 2.09f, 1.1f, 4.13f, 3.31f, 2.8f, 5.83f, 5.09f, 6.1f, 2.86f, 0.24f, 5.76f, 5.19f, 0.72f
};
//! 1 charge
static auto const c_sampleCharges1 = ChargesVector::fromVector(c_sampleChargesFull.begin(), c_sampleChargesFull.begin() + 1);
//! 2 charges
static auto const c_sampleCharges2 = ChargesVector::fromVector(c_sampleChargesFull.begin() + 1, c_sampleChargesFull.begin() + 3);
//! 13 charges
static auto const c_sampleCharges13 = ChargesVector::fromVector(c_sampleChargesFull.begin() + 3, c_sampleChargesFull.begin() + 16);

//! Random coordinate vectors
static CoordinatesVector const c_sampleCoordinatesFull
{
    {
        5.59f, 1.37f, 0.95f
    }, {
        16.0f, 1.02f, 0.22f      // 2 box lengths in x
    }, {
        0.034f, 1.65f, 0.22f
    }, {
        0.33f, 0.92f, 1.56f
    }, {
        1.16f, 0.75f, 0.39f
    }, {
        0.5f, 1.63f, 1.14f
    }, {
        16.0001f, 1.52f, 1.19f  // > 2 box lengths in x
    }, {
        1.43f, 1.1f, 4.1f       // > 2 box lengths in z
    }, {
        -1.08f, 1.19f, 0.08f    // negative x
    }, {
        1.6f, 0.93f, 0.53f
    }, {
        1.32f, -1.48f, 0.16f    // negative y
    }, {
        0.87f, 0.0f, 0.33f
    }, {
        0.95f, 7.7f, -0.48f     // > 2 box lengths in y, negative z
    }, {
        1.23f, 0.91f, 0.68f
    }, {
        0.19f, 1.45f, 0.94f
    }, {
        1.28f, 0.46f, 0.38f
    }, {
        1.21f, 0.23f, 1.0f
    }
};
//! 1 coordinate vector
static CoordinatesVector const c_sampleCoordinates1(c_sampleCoordinatesFull.begin(), c_sampleCoordinatesFull.begin() + 1);
//! 2 coordinate vectors
static CoordinatesVector const c_sampleCoordinates2(c_sampleCoordinatesFull.begin() + 1, c_sampleCoordinatesFull.begin() + 3);
//! 13 coordinate vectors
static CoordinatesVector const c_sampleCoordinates13(c_sampleCoordinatesFull.begin() + 3, c_sampleCoordinatesFull.begin() + 16);

//! moved out from instantiantions for readability
auto c_inputBoxes     = ::testing::ValuesIn(c_sampleBoxes);
//! moved out from instantiantions for readability
auto c_inputPmeOrders = ::testing::Range(3, 5 + 1);
//! moved out from instantiantions for readability
auto c_inputGridSizes = ::testing::ValuesIn(c_sampleGridSizes);

/*! \brief Instantiation of the test with valid input and 1 atom */
INSTANTIATE_TEST_CASE_P(SaneInput1, PmeSplineAndSpreadTest, ::testing::Combine(c_inputBoxes, c_inputPmeOrders, c_inputGridSizes,
                                                                                   ::testing::Values(c_sampleCoordinates1),
                                                                                   ::testing::Values(c_sampleCharges1)
                                                                               ));

/*! \brief Instantiation of the test with valid input and 2 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput2, PmeSplineAndSpreadTest, ::testing::Combine(c_inputBoxes, c_inputPmeOrders, c_inputGridSizes,
                                                                                   ::testing::Values(c_sampleCoordinates2),
                                                                                   ::testing::Values(c_sampleCharges2)
                                                                               ));
/*! \brief Instantiation of the test with valid input and 13 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput13, PmeSplineAndSpreadTest, ::testing::Combine(c_inputBoxes, c_inputPmeOrders, c_inputGridSizes,
                                                                                    ::testing::Values(c_sampleCoordinates13),
                                                                                    ::testing::Values(c_sampleCharges13)
                                                                                ));
}
}
}
