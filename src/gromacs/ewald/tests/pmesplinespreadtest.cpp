/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

#include "config.h"

#include <string>

#include <gmock/gmock.h>

#include "gromacs/mdrunutility/mdmodules.h"
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
    private:
        //! Environment for getting the t_inputrec structure easily
        MDModules mdModules_;
    public:
        //! Default constructor
        PmeSplineAndSpreadTest() = default;
        //! The test
        void RunTest()
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
            t_inputrec *inputRec  = mdModules_.inputrec();
            inputRec->nkx         = gridSize[XX];
            inputRec->nky         = gridSize[YY];
            inputRec->nkz         = gridSize[ZZ];
            inputRec->pme_order   = pmeOrder;
            inputRec->coulombtype = eelPME;
            inputRec->epsilon_r   = 1.0;

            TestReferenceData                                      refData;

            const std::map<PmeCodePath, std::string>               modesToTest = {{PmeCodePath::CPU, "CPU"},
#if GMX_GPU == GMX_GPU_CUDA
                                                                                  {PmeCodePath::CUDA, "CUDA"},
#endif
            };

            const std::map<PmeSplineAndSpreadOptions, std::string> optionsToTest = {{PmeSplineAndSpreadOptions::SplineAndSpreadUnified, "spline computation and charge spreading (fused)"},
                                                                                    {PmeSplineAndSpreadOptions::SplineOnly, "spline computation"},
                                                                                    {PmeSplineAndSpreadOptions::SpreadOnly, "charge spreading"}};
            for (const auto &mode : modesToTest)
            {
                for (const auto &option : optionsToTest)
                {
                    /* Describing the test uniquely in case it fails */

                    SCOPED_TRACE(formatString("Testing %s with %s for PME grid size %d %d %d"
                                              ", order %d, %zu atoms",
                                              option.second.c_str(), mode.second.c_str(),
                                              gridSize[XX], gridSize[YY], gridSize[ZZ],
                                              pmeOrder,
                                              atomCount));

                    /* Running the test */

                    const bool supportedInput = PmeSupportsInputForMode(inputRec, mode.first);
                    if (!supportedInput)
                    {
                        EXPECT_THROW(PmeInitWithAtoms(inputRec, mode.first, coordinates, charges, box), NotImplementedError);
                    }
                    else
                    {
                        PmeSafePointer pmeSafe(PmeInitWithAtoms(inputRec, mode.first, coordinates, charges, box));

                        const bool     computeSplines = (option.first == PmeSplineAndSpreadOptions::SplineOnly) || (option.first == PmeSplineAndSpreadOptions::SplineAndSpreadUnified);
                        const bool     spreadCharges  = (option.first == PmeSplineAndSpreadOptions::SpreadOnly) || (option.first == PmeSplineAndSpreadOptions::SplineAndSpreadUnified);

                        if (!computeSplines)
                        {
                            // Here we should set up the results of the spline computation so that the spread can run.
                            // What is lazy and works is running the separate spline so that it will set it up for us:
                            PmePerformSplineAndSpread(pmeSafe, mode.first, true, false);
                            // We know that it is tested in another iteration.
                            // TODO: Clean alternative: read and set the reference gridline indices, spline params
                            // - with currently missing reference checker readSequence methods.
                        }

                        PmePerformSplineAndSpread(pmeSafe, mode.first, computeSplines, spreadCharges);

                        /* Outputs correctness check */
                        /* All tolerances were picked empirically for single precision on CPU */

                        TestReferenceChecker checker(refData.rootChecker());

                        if (computeSplines)
                        {
                            SplineParamsVector    splineValues, splineDerivatives;
                            GridLineIndicesVector gridLineIndices;
                            PmeFetchOutputsSpline(pmeSafe, mode.first, splineValues, splineDerivatives, gridLineIndices);

                            /* Spline data - values and derivatives */
                            const auto            ulpToleranceSplines = 40;
                            checker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpToleranceSplines));
                            checker.checkSequence(splineValues.begin(), splineValues.end(), "Spline values");
                            checker.checkSequence(splineDerivatives.begin(), splineDerivatives.end(), "Spline derivatives");

                            /* Particle gridline indices */
                            checker.checkSequence(gridLineIndices.begin(), gridLineIndices.end(), "Gridline indices");
                        }

                        if (spreadCharges)
                        {
                            SparseGridValues nonZeroGridValues;
                            PmeFetchOutputsSpread(pmeSafe, mode.first, nonZeroGridValues);

                            /* The wrapped grid */
                            TestReferenceChecker gridValuesChecker(checker.checkCompound("NonZeroGridValues", "RealSpaceGrid"));
                            const auto           ulpTolerance = 48;
                            gridValuesChecker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpTolerance));
                            for (const auto &point : nonZeroGridValues)
                            {
                                std::string valueId = formatString("Cell %d %d %d", point.first[XX], point.first[YY], point.first[ZZ]);
                                gridValuesChecker.checkReal(point.second, valueId.c_str());
                            }
                        }
                    }
                }
            }
        }
};


/*! \brief Test for PME B-spline moduli computation */
TEST_P(PmeSplineAndSpreadTest, ReproducesOutputs)
{
    EXPECT_NO_THROW(RunTest());
}

/* Valid input instances */

//! A couple of valid inputs for boxes.
static std::vector<Matrix3x3> const sampleBoxes
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
                   0.0f, 2.1f, 0.0f,
                   6.0f, 5.0f, 3.2f
               }},
};

//! A couple of valid inputs for grid sizes.
static std::vector<IVec> const sampleGridSizes
{
    IVec {
        16, 12, 14
    },
    IVec {
        9, 7, 11
    }
};

//! Random charges
static ChargesVector const sampleChargesFull
{
    4.95f, 3.11f, 3.97f, 1.08f, 2.09f, 1.1f, 4.13f, 3.31f, 2.8f, 5.83f, 5.09f, 6.1f, 2.86f, 0.24f, 5.76f, 5.19f, 0.72f
};
//! 1 charge
static ChargesVector const sampleCharges1(sampleChargesFull.begin(), sampleChargesFull.begin() + 1);
//! 2 charges
static ChargesVector const sampleCharges2(sampleChargesFull.begin() + 1, sampleChargesFull.begin() + 3);
//! 13 charges
static ChargesVector const sampleCharges13(sampleChargesFull.begin() + 3, sampleChargesFull.begin() + 16);

//! Random coordinate vectors
static CoordinatesVector const sampleCoordinatesFull
{
    {
        1.59f, 1.37f, 0.95f
    }, {
        16.0f, 1.02f, 0.22f  // 2 box lengths in x
    }, {
        0.03f, 1.6f, 0.22f
    }, {
        0.33f, 0.92f, 1.56f
    }, {
        1.16f, 0.75f, 0.39f
    }, {
        0.49f, 1.6f, 1.12f
    }, {
        16.0001f, 1.52f, 1.19f  // >2 box lengths in x
    }, {
        1.43f, 1.1f, 0.08f
    }, {
        -1.08f, 1.19f, 0.08f   // negative x
    }, {
        1.6f, 0.93f, 0.53f
    }, {
        1.32f, 1.48f, 0.16f
    }, {
        0.87f, 0.0f, 0.33f
    }, {
        0.95f, 0.13f, 0.48f
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
static CoordinatesVector const sampleCoordinates1(sampleCoordinatesFull.begin(), sampleCoordinatesFull.begin() + 1);
//! 2 coordinate vectors
static CoordinatesVector const sampleCoordinates2(sampleCoordinatesFull.begin() + 1, sampleCoordinatesFull.begin() + 3);
//! 13 coordinate vectors
static CoordinatesVector const sampleCoordinates13(sampleCoordinatesFull.begin() + 3, sampleCoordinatesFull.begin() + 16);

//! moved out from instantiantions for readability
auto inputBoxes     = ::testing::ValuesIn(sampleBoxes);
//! moved out from instantiantions for readability
auto inputPmeOrders = ::testing::Range(4, 5 + 1);
//! moved out from instantiantions for readability
auto inputGridSizes = ::testing::ValuesIn(sampleGridSizes);

/*! \brief Instantiation of the test with valid input and 1 atom */
INSTANTIATE_TEST_CASE_P(SaneInput1, PmeSplineAndSpreadTest, ::testing::Combine(inputBoxes, inputPmeOrders, inputGridSizes,
                                                                                   ::testing::Values(sampleCoordinates1),
                                                                                   ::testing::Values(sampleCharges1)
                                                                               ));

/*! \brief Instantiation of the test with valid input and 2 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput2, PmeSplineAndSpreadTest, ::testing::Combine(inputBoxes, inputPmeOrders, inputGridSizes,
                                                                                   ::testing::Values(sampleCoordinates2),
                                                                                   ::testing::Values(sampleCharges2)
                                                                               ));
/*! \brief Instantiation of the test with valid input and 13 atoms */
INSTANTIATE_TEST_CASE_P(SaneInput13, PmeSplineAndSpreadTest, ::testing::Combine(inputBoxes, inputPmeOrders, inputGridSizes,
                                                                                    ::testing::Values(sampleCoordinates13),
                                                                                    ::testing::Values(sampleCharges13)
                                                                                ));
}
}
}
