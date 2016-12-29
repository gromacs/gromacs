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
 * Implements PME force gathering tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

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
/*! \brief Convenience typedef of the test input parameters - unit cell box, PME interpolation order, grid dimensions,
 *  grid values, overwriting/reducing the input forces, gridline indices, spline theta values, spline dtheta values, atom charges,
 */
typedef std::tuple<Matrix3x3, int, IVec, SparseGridValues, PmeGatherInputHandling, GridLineIndicesVector, SplineParamsVector,
                   SplineParamsVector, ChargesVector> GatherInputParameters;

//! Test fixture
class PmeGatherTest : public ::testing::TestWithParam<GatherInputParameters>
{
    private:
        //! Environment for getting the t_inputrec structure easily
        MDModules mdModules_;

    public:
        //! Default constructor
        PmeGatherTest() = default;

        //! The test
        void RunTest()
        {
            /* Getting the input */
            Matrix3x3                 box;
            int                       pmeOrder;
            IVec                      gridSize;
            ChargesVector             charges;
            GridLineIndicesVector     gridLineIndices;
            SplineParamsVector        splineValues;
            SplineParamsVector        splineDerivatives;
            SparseGridValues          nonZeroGridValues;
            PmeGatherInputHandling    inputForceTreatment;
            std::tie(box, pmeOrder, gridSize, nonZeroGridValues, inputForceTreatment, gridLineIndices, splineValues, splineDerivatives, charges) = GetParam();
            const size_t              atomCount = charges.size();
            CoordinatesVector         coordinatesDummy(atomCount, RVec {1e6, 1e7, -1e8});
            /* The coordinates are intentionaly bogus - only the size matters; the gridline indices are fed directly as inputs */

            /* Storing the input where it's needed, running the test */
            t_inputrec *inputRec  = mdModules_.inputrec();
            inputRec->nkx         = gridSize[XX];
            inputRec->nky         = gridSize[YY];
            inputRec->nkz         = gridSize[ZZ];
            inputRec->pme_order   = pmeOrder;
            inputRec->coulombtype = eelPME;

            TestReferenceData refData;
            const std::map<PmeCodePath, std::string> modesToTest = {{PmeCodePath::CPU, "CPU"}};
            for (const auto &mode : modesToTest)
            {
                /* Describing the test uniquely */
                SCOPED_TRACE(formatString("Testing force gathering with %s for PME grid size %d %d %d"
                                          ", order %d, %zu atoms",
                                          mode.second.c_str(),
                                          gridSize[XX], gridSize[YY], gridSize[ZZ],
                                          pmeOrder,
                                          atomCount));

                /* Running the test */
                PmeSafePointer  pmeSafe = PmeInitWithAtoms(inputRec, coordinatesDummy, charges, box);

                PmeSetRealGrid(pmeSafe, mode.first, nonZeroGridValues);
                PmeSetGridLineIndices(pmeSafe, mode.first, gridLineIndices);
                PmeSetSplineData(pmeSafe, mode.first, splineValues, PmeSplineDataType::Values);
                PmeSetSplineData(pmeSafe, mode.first, splineDerivatives, PmeSplineDataType::Derivatives);

                const ForcesVector inputForcesFull {{0.02, 0.87, 0.95}, {0.66, 0.67, 0.38}, {0.45, 0.04, 0.94}, {0.54, 0.76, 0.58}, {0.83, 0.31, 0.73}, {0.71, 0.06, 0.35}, {0.32, 0.35, 0.61}, {0.27, 0.98, 0.83}, {0.11, 0.3, 0.42}, {0.95, 0.69, 0.58}, {0.29, 0.1, 0.68}, {0.94, 0.62, 0.51}, {0.47, 0.04, 0.47}, {0.34, 0.71, 0.52}};
                GMX_RELEASE_ASSERT(inputForcesFull.size() >= atomCount, "Bad input forces size");
                ForcesVector       forces(inputForcesFull.begin(), inputForcesFull.begin() + atomCount);
                PmePerformGather(pmeSafe, mode.first, inputForceTreatment, forces);

                /* Check the output forces correctness */
                TestReferenceChecker checker(refData.rootChecker());
                const auto           ulpTolerance = 4;
                checker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpTolerance));
                checker.checkSequence(forces.begin(), forces.end(), "Forces");
            }
        }
};

// TODO fix the same comment on the spreading test
/*! \brief Test for PME force gathering */
TEST_P(PmeGatherTest, ReproducesOutputs)
{
    EXPECT_NO_THROW(RunTest());
}

// copypasted from the spreading tests. TODO: reuse

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
                   0.0f, 3.1f, 0.0f,
                   6.0f, 4.0f, 3.2f
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

//! moved out from instantiations for readability
auto inputBoxes     = ::testing::ValuesIn(sampleBoxes);
//! moved out from instantiations for readability
auto inputGridSizes = ::testing::ValuesIn(sampleGridSizes);

// COPYPASTA END

//! All the input atom gridline indices
static GridLineIndicesVector const sampleGridLineIndicesFull
{
    IVec {
        4, 2, 6
    },
    IVec {
        1, 4, 10
    },
    IVec {
        0, 6, 6
    },
    IVec {
        0, 1, 4
    },
    IVec {
        6, 3, 0
    },
    IVec {
        7, 2, 2
    },
    IVec {
        8, 3, 1
    },
    IVec {
        4, 0, 3
    },
    IVec {
        0, 0, 0
    },
    IVec {
        8, 5, 8
    },
    IVec {
        4, 4, 2
    },
    IVec {
        7, 1, 7
    },
    IVec {
        8, 5, 5
    },
    IVec {
        2, 6, 5
    },
    IVec {
        1, 6, 2
    },
    IVec {
        7, 1, 8
    },
    IVec {
        3, 5, 1
    },
};
//! 1 atom
static GridLineIndicesVector const sampleGridLineIndices1(sampleGridLineIndicesFull.begin(), sampleGridLineIndicesFull.begin() + 1);
//! 2 atoms
static GridLineIndicesVector const sampleGridLineIndices2(sampleGridLineIndicesFull.begin() + 1, sampleGridLineIndicesFull.begin() + 3);
//! 13 atoms
static GridLineIndicesVector const sampleGridLineIndices13(sampleGridLineIndicesFull.begin() + 3, sampleGridLineIndicesFull.begin() + 16);

// Spline values/derivatives below are also generated randomly, so they are bogus, but that should not affect the reproducibility, which we're after
//! Helper constant for spline data striding, PME order 4
const auto splineDataStride4 = DIM * 4;
//! Helper constant for spline data striding, PME order 5
const auto splineDataStride5 = DIM * 5;

//! A lot of bogus input spline values - should have at list (max PME order = 5) * (DIM = 3) * (total unique atom number in all test cases = 16) values
static SplineParamsVector const sampleSplineValuesFull
{
    0.12, 0.81, 0.29, 0.22, 0.13, 0.19, 0.12, 0.8, 0.44, 0.38, 0.32, 0.36, 0.27, 0.11, 0.17, 0.94, 0.07, 0.9, 0.98, 0.96, 0.07, 0.94, 0.77, 0.24, 0.84, 0.16, 0.77, 0.57, 0.52, 0.27, 0.39, 0.45, 0.6, 0.59, 0.44, 0.91, 0.97, 0.43, 0.24, 0.52, 0.73, 0.55, 0.99, 0.39, 0.97, 0.35, 0.1, 0.68, 0.19, 0.1, 0.77, 0.2, 0.43, 0.69, 0.76, 0.32, 0.31, 0.94, 0.53, 0.6, 0.93, 0.57, 0.94, 0.88, 0.75, 0.77, 0.91, 0.72, 0.07, 0.78, 0.09, 0.02, 0.48, 0.97, 0.89, 0.39, 0.48, 0.19, 0.02, 0.92, 0.8, 0.41, 0.53, 0.32, 0.38, 0.58, 0.36, 0.46, 0.92, 0.91, 0.01, 0.86, 0.54, 0.86, 0.94, 0.37, 0.35, 0.81, 0.89, 0.48,
    0.34, 0.18, 0.11, 0.02, 0.87, 0.95, 0.66, 0.67, 0.38, 0.45, 0.04, 0.94, 0.54, 0.76, 0.58, 0.83, 0.31, 0.73, 0.71, 0.06, 0.35, 0.32, 0.35, 0.61, 0.27, 0.98, 0.83, 0.11, 0.3, 0.42, 0.95, 0.69, 0.58, 0.29, 0.1, 0.68, 0.94, 0.62, 0.51, 0.47, 0.04, 0.47, 0.34, 0.71, 0.52, 0.19, 0.69, 0.5, 0.59, 0.05, 0.74, 0.11, 0.4, 0.81, 0.24, 0.53, 0.71, 0.07, 0.17, 0.41, 0.23, 0.78, 0.27, 0.1, 0.71, 0.36, 0.67, 0.6, 0.94, 0.69, 0.19, 0.58, 0.68, 0.5, 0.62, 0.38, 0.29, 0.44, 0.04, 0.89, 0, 0.76, 0.22, 0.16, 0.08, 0.62, 0.51, 0.62, 0.83, 0.72, 0.96, 0.99, 0.4, 0.79, 0.83, 0.21, 0.43, 0.32, 0.44, 0.72,
    0.21, 0.4, 0.93, 0.07, 0.11, 0.41, 0.24, 0.04, 0.36, 0.15, 0.92, 0.08, 0.99, 0.35, 0.42, 0.7, 0.17, 0.39, 0.69, 0, 0.86, 0.89, 0.59, 0.81, 0.77, 0.15, 0.89, 0.17, 0.76, 0.67, 0.58, 0.78, 0.26, 0.19, 0.69, 0.18, 0.46, 0.6, 0.69, 0.23, 0.34, 0.3, 0.64, 0.34, 0.6, 0.99, 0.69, 0.57, 0.75, 0.07, 0.36, 0.75, 0.81, 0.8, 0.42, 0.09, 0.94, 0.66, 0.35, 0.67, 0.34, 0.66, 0.02, 0.47, 0.78, 0.21, 0.02, 0.18, 0.42, 0.2, 0.46, 0.34, 0.4, 0.46, 0.96, 0.86, 0.25, 0.25, 0.22, 0.37, 0.59, 0.19, 0.45, 0.61, 0.04, 0.71, 0.77, 0.51, 0.77, 0.15, 0.78, 0.36, 0.62, 0.24, 0.86, 0.2, 0.77, 0.08, 0.09, 0.3,
    0.0, 0.6, 0.99, 0.69,
};
//! 1 atom, order 4
static SplineParamsVector const sampleSplineValues4_1(sampleSplineValuesFull.begin(),
                                                      sampleSplineValuesFull.begin() + splineDataStride4 * 1);
//! 2 atoms, order 4
static SplineParamsVector const sampleSplineValues4_2(sampleSplineValuesFull.begin() + splineDataStride4 * 1,
                                                      sampleSplineValuesFull.begin() + splineDataStride4 * 3);
//! 13 atoms, order 4
static SplineParamsVector const sampleSplineValues4_13(sampleSplineValuesFull.begin() + splineDataStride4 * 3,
                                                       sampleSplineValuesFull.begin() + splineDataStride4 * 16);
//! 1 atom, order 5
static SplineParamsVector const sampleSplineValues5_1(sampleSplineValuesFull.begin(),
                                                      sampleSplineValuesFull.begin() + splineDataStride5 * 1);
//! 2 atoms, order 5
static SplineParamsVector const sampleSplineValues5_2(sampleSplineValuesFull.begin() + splineDataStride5 * 1,
                                                      sampleSplineValuesFull.begin() + splineDataStride5 * 3);
//! 13 atoms, order 5
static SplineParamsVector const sampleSplineValues5_13(sampleSplineValuesFull.begin() + splineDataStride5 * 3,
                                                       sampleSplineValuesFull.begin() + splineDataStride5 * 16);

//! A lot of bogus input spline derivatives - should have at list (max PME order = 5) * (DIM = 3) * (total unique atom number in all test cases = 16) values
static SplineParamsVector const sampleSplineDerivativesFull
{
    0.82, 0.88, 0.83, 0.11, 0.93, 0.32, 0.71, 0.37, 0.69, 0.88, 0.11, 0.38, 0.25, 0.5, 0.36, 0.81, 0.78, 0.31, 0.66, 0.32, 0.27, 0.35, 0.53, 0.83, 0.08, 0.08, 0.94, 0.71, 0.65, 0.24, 0.13, 0.01, 0.33, 0.65, 0.24, 0.53, 0.45, 0.84, 0.33, 0.97, 0.31, 0.7, 0.03, 0.31, 0.41, 0.76, 0.12, 0.3, 0.57, 0.65, 0.87, 0.99, 0.42, 0.97, 0.32, 0.39, 0.73, 0.23, 0.03, 0.67, 0.97, 0.57, 0.42, 0.38, 0.54, 0.17, 0.53, 0.54, 0.18, 0.8, 0.76, 0.13, 0.29, 0.83, 0.77, 0.56, 0.4, 0.87, 0.36, 0.18, 0.59, 0.04, 0.05, 0.61, 0.26, 0.91, 0.62, 0.16, 0.89, 0.23, 0.26, 0.59, 0.33, 0.2, 0.49, 0.41, 0.25, 0.4, 0.16, 0.83,
    0.44, 0.82, 0.21, 0.95, 0.14, 0.8, 0.37, 0.31, 0.41, 0.53, 0.15, 0.85, 0.78, 0.17, 0.92, 0.03, 0.13, 0.2, 0.03, 0.33, 0.87, 0.38, 0, 0.08, 0.79, 0.36, 0.53, 0.05, 0.07, 0.94, 0.23, 0.85, 0.13, 0.27, 0.23, 0.22, 0.26, 0.38, 0.15, 0.48, 0.18, 0.33, 0.23, 0.62, 0.1, 0.36, 0.99, 0.07, 0.02, 0.04, 0.09, 0.29, 0.52, 0.29, 0.83, 0.97, 0.61, 0.81, 0.49, 0.56, 0.08, 0.09, 0.03, 0.65, 0.46, 0.1, 0.06, 0.06, 0.39, 0.29, 0.04, 0.03, 0.1, 0.83, 0.94, 0.59, 0.97, 0.82, 0.2, 0.66, 0.23, 0.11, 0.03, 0.16, 0.27, 0.53, 0.94, 0.46, 0.43, 0.29, 0.97, 0.64, 0.46, 0.37, 0.43, 0.48, 0.37, 0.93, 0.5, 0.2,
    0.92, 0.09, 0.74, 0.55, 0.44, 0.05, 0.13, 0.17, 0.79, 0.44, 0.11, 0.6, 0.64, 0.05, 0.96, 0.3, 0.45, 0.47, 0.42, 0.74, 0.91, 0.06, 0.89, 0.24, 0.26, 0.68, 0.4, 0.88, 0.5, 0.65, 0.48, 0.15, 0, 0.41, 0.67, 0.4, 0.31, 0.73, 0.77, 0.36, 0.26, 0.74, 0.46, 0.56, 0.78, 0.92, 0.32, 0.9, 0.06, 0.55, 0.6, 0.13, 0.38, 0.93, 0.5, 0.92, 0.96, 0.82, 0, 0.04, 0.9, 0.55, 0.97, 1, 0.23, 0.46, 0.52, 0.49, 0, 0.32, 0.16, 0.4, 0.62, 0.36, 0.03, 0.63, 0.16, 0.58, 0.97, 0.03, 0.44, 0.07, 0.22, 0.75, 0.32, 0.61, 0.94, 0.33, 0.7, 0.57, 0.5, 0.84, 0.7, 0.47, 0.18, 0.09, 0.25, 0.77, 0.94, 0.85,
    0.09, 0.83, 0.02, 0.91,
};
//! 1 atom, order 4
static SplineParamsVector const sampleSplineDerivatives4_1(sampleSplineDerivativesFull.begin(),
                                                           sampleSplineDerivativesFull.begin() + splineDataStride4 * 1);
//! 2 atoms, order 4
static SplineParamsVector const sampleSplineDerivatives4_2(sampleSplineDerivativesFull.begin() + splineDataStride4 * 1,
                                                           sampleSplineDerivativesFull.begin() + splineDataStride4 * 3);
//! 13 atoms, order 4
static SplineParamsVector const sampleSplineDerivatives4_13(sampleSplineDerivativesFull.begin() + splineDataStride4 * 3,
                                                            sampleSplineDerivativesFull.begin() + splineDataStride4 * 16);
//! 1 atom, order 5
static SplineParamsVector const sampleSplineDerivatives5_1(sampleSplineDerivativesFull.begin(),
                                                           sampleSplineDerivativesFull.begin() + splineDataStride5 * 1);
//! 2 atoms, order 5
static SplineParamsVector const sampleSplineDerivatives5_2(sampleSplineDerivativesFull.begin() + splineDataStride5 * 1,
                                                           sampleSplineDerivativesFull.begin() + splineDataStride5 * 3);
//! 13 atoms, order 5
static SplineParamsVector const sampleSplineDerivatives5_13(sampleSplineDerivativesFull.begin() + splineDataStride5 * 3,
                                                            sampleSplineDerivativesFull.begin() + splineDataStride5 * 16);

//! 2 sample grids - only non-zero values have to be listed
static std::vector<SparseGridValues> const sampleGrids
{
    SparseGridValues {{
                          IVec {
                              0, 0, 0
                          }, 3.5f
                      }, {
                          IVec {
                              7, 0, 0
                          }, -2.5f
                      }, {
                          IVec {
                              3, 5, 7
                          }, -0.006f
                      }, {
                          IVec {
                              3, 1, 2
                          }, 0.6f
                      },  {
                          IVec {
                              6, 2, 4
                          }, 300.1f
                      }, },
    SparseGridValues {{
                          IVec {
                              0, 4, 0
                          }, 6.f
                      }, {
                          IVec {
                              4, 2, 7
                          }, 13.76f
                      }, {
                          IVec {
                              0, 6, 7
                          }, 3.6f
                      }, {
                          IVec {
                              2, 5, 10
                          }, 3.6f
                      }, }
};

//! Moved out from instantiations for readability
auto inputGrids = ::testing::ValuesIn(sampleGrids);
//! Moved out from instantiations for readability
auto instPmeOrder4 = ::testing::Values(4);
//! Moved out from instantiations for readability
auto instPmeOrder5 = ::testing::Values(5);

//! Moved out from instantiations for readability
auto inputForceTreatment = ::testing::Values(PmeGatherInputHandling::Overwrite, PmeGatherInputHandling::ReduceWith);

//! Instantiation of the PME gathering test with 1 atom, order of 4
INSTANTIATE_TEST_CASE_P(SaneInput1_4, PmeGatherTest, ::testing::Combine(inputBoxes, instPmeOrder4, inputGridSizes, inputGrids, inputForceTreatment,
                                                                            ::testing::Values(sampleGridLineIndices1),
                                                                            ::testing::Values(sampleSplineValues4_1),
                                                                            ::testing::Values(sampleSplineDerivatives4_1),
                                                                            ::testing::Values(sampleCharges1)
                                                                        ));
//! Instantiation of the PME gathering test with 2 atoms, order of 4
INSTANTIATE_TEST_CASE_P(SaneInput2_4, PmeGatherTest, ::testing::Combine(inputBoxes, instPmeOrder4, inputGridSizes, inputGrids, inputForceTreatment,
                                                                            ::testing::Values(sampleGridLineIndices2),
                                                                            ::testing::Values(sampleSplineValues4_2),
                                                                            ::testing::Values(sampleSplineDerivatives4_2),
                                                                            ::testing::Values(sampleCharges2)
                                                                        ));
//! Instantiation of the PME gathering test with 13 atoms, order of 4
INSTANTIATE_TEST_CASE_P(SaneInput13_4, PmeGatherTest, ::testing::Combine(inputBoxes, instPmeOrder4, inputGridSizes, inputGrids, inputForceTreatment,
                                                                             ::testing::Values(sampleGridLineIndices13),
                                                                             ::testing::Values(sampleSplineValues4_13),
                                                                             ::testing::Values(sampleSplineDerivatives4_13),
                                                                             ::testing::Values(sampleCharges13)
                                                                         ));
//! Instantiation of the PME gathering test with 1 atom, order of 5
INSTANTIATE_TEST_CASE_P(SaneInput1_5, PmeGatherTest, ::testing::Combine(inputBoxes, instPmeOrder5, inputGridSizes, inputGrids, inputForceTreatment,
                                                                            ::testing::Values(sampleGridLineIndices1),
                                                                            ::testing::Values(sampleSplineValues5_1),
                                                                            ::testing::Values(sampleSplineDerivatives5_1),
                                                                            ::testing::Values(sampleCharges1)
                                                                        ));
//! Instantiation of the PME gathering test with 2 atoms, order of 5
INSTANTIATE_TEST_CASE_P(SaneInput2_5, PmeGatherTest, ::testing::Combine(inputBoxes, instPmeOrder5, inputGridSizes, inputGrids, inputForceTreatment,
                                                                            ::testing::Values(sampleGridLineIndices2),
                                                                            ::testing::Values(sampleSplineValues5_2),
                                                                            ::testing::Values(sampleSplineDerivatives5_2),
                                                                            ::testing::Values(sampleCharges2)
                                                                        ));
//! Instantiation of the PME gathering test with 13 atoms, order of 5
INSTANTIATE_TEST_CASE_P(SaneInput13_5, PmeGatherTest, ::testing::Combine(inputBoxes, instPmeOrder5, inputGridSizes, inputGrids, inputForceTreatment,
                                                                             ::testing::Values(sampleGridLineIndices13),
                                                                             ::testing::Values(sampleSplineValues5_13),
                                                                             ::testing::Values(sampleSplineDerivatives5_13),
                                                                             ::testing::Values(sampleCharges13)
                                                                         ));

}
}
}
