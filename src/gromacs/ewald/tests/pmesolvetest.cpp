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
 * Implements PME solving tests.
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
/*! \brief Convenience typedef of the test input parameters - unit cell box, complex grid dimensions, complex grid values,
 * electrostatic constant epsilon_r, Ewald splitting parameters ewaldcoeff_q and ewaldcoeff_lj, solver type
 * Output: transformed local grid (Fourier space); optionally reciprocal energy and virial matrix.
 * TODO:
 * Implement and test Lorentz-Berthelot
 */
typedef std::tuple<Matrix3x3, IVec, SparseComplexGridValuesInput, double, double, double, PmeSolveAlgorithm> SolveInputParameters;

//! Test fixture
class PmeSolveTest : public ::testing::TestWithParam<SolveInputParameters>
{
    public:
        //! Default constructor
        PmeSolveTest() = default;

        //! The test
        void runTest()
        {
            /* Getting the input */
            Matrix3x3                      box;
            IVec                           gridSize;
            SparseComplexGridValuesInput   nonZeroGridValues;
            double                         epsilon_r;
            double                         ewaldCoeff_q;
            double                         ewaldCoeff_lj;
            PmeSolveAlgorithm              method;
            std::tie(box, gridSize, nonZeroGridValues, epsilon_r, ewaldCoeff_q, ewaldCoeff_lj, method) = GetParam();

            /* Storing the input where it's needed, running the test */
            t_inputrec inputRec;
            inputRec.nkx         = gridSize[XX];
            inputRec.nky         = gridSize[YY];
            inputRec.nkz         = gridSize[ZZ];
            inputRec.pme_order   = 4;
            inputRec.coulombtype = eelPME;
            inputRec.epsilon_r   = epsilon_r;
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb:
                    break;

                case PmeSolveAlgorithm::LennardJones:
                    inputRec.vdwtype = evdwPME;
                    break;

                default:
                    GMX_THROW(InternalError("Unknown PME solver"));
            }

            TestReferenceData                     refData;
            const std::map<CodePath, std::string> modesToTest = {{CodePath::CPU, "CPU"},
                                                                 {CodePath::CUDA, "CUDA"}};
            for (const auto &mode : modesToTest)
            {
                const bool supportedInput = pmeSupportsInputForMode(&inputRec, mode.first);
                if (!supportedInput)
                {
                    /* Testing the failure for the unsupported input */
                    EXPECT_THROW(pmeInitEmpty(&inputRec, mode.first, nullptr, box, ewaldCoeff_q, ewaldCoeff_lj), NotImplementedError);
                    continue;
                }

                std::map<GridOrdering, std::string> gridOrderingsToTest = {{GridOrdering::YZX, "YZX"}};
                if (mode.first == CodePath::CUDA)
                {
                    gridOrderingsToTest[GridOrdering::XYZ] = "XYZ";
                }
                const auto contextsToTest = pmeEnv->getHardwareContexts(mode.first);
                for (const auto &gridOrdering : gridOrderingsToTest)
                {
                    for (const auto &context : contextsToTest)
                    {
                        for (bool computeEnergyAndVirial : {false, true})
                        {
                            /* Describing the test*/
                            SCOPED_TRACE(formatString("Testing solving (%s, %s, %s energy/virial) with %s %sfor PME grid size %d %d %d, Ewald coefficients %g %g",
                                                      (method == PmeSolveAlgorithm::LennardJones) ? "Lennard-Jones" : "Coulomb",
                                                      gridOrdering.second.c_str(),
                                                      computeEnergyAndVirial ? "with" : "without",
                                                      mode.second.c_str(),
                                                      context.getDescription().c_str(),
                                                      gridSize[XX], gridSize[YY], gridSize[ZZ],
                                                      ewaldCoeff_q, ewaldCoeff_lj
                                                      ));

                            /* Running the test */
                            PmeSafePointer pmeSafe = pmeInitEmpty(&inputRec, mode.first, context.getDeviceInfo(), box, ewaldCoeff_q, ewaldCoeff_lj);
                            pmeSetComplexGrid(pmeSafe.get(), mode.first, gridOrdering.first, nonZeroGridValues);
                            const real     cellVolume = box[0] * box[4] * box[8];
                            //FIXME - this is box[XX][XX] * box[YY][YY] * box[ZZ][ZZ], should be stored in the PME structure
                            pmePerformSolve(pmeSafe.get(), mode.first, method, cellVolume, gridOrdering.first, computeEnergyAndVirial);
                            pmeFinalizeTest(pmeSafe.get(), mode.first);

                            /* Check the outputs */
                            TestReferenceChecker          checker(refData.rootChecker());

                            SparseComplexGridValuesOutput nonZeroGridValuesOutput = pmeGetComplexGrid(pmeSafe.get(), mode.first, gridOrdering.first);
                            /* Transformed grid */
                            TestReferenceChecker          gridValuesChecker(checker.checkCompound("NonZeroGridValues", "ComplexSpaceGrid"));
                            const auto                    ulpToleranceGrid = 40;
                            gridValuesChecker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpToleranceGrid));
                            for (const auto &point : nonZeroGridValuesOutput)
                            {
                                // we want an additional safeguard for denormal numbers as they cause an exception in string conversion;
                                // however, using GMX_REAL_MIN causes an "unused item warning" for single precision builds
                                if (fabs(point.second.re) >= GMX_FLOAT_MIN)
                                {
                                    gridValuesChecker.checkReal(point.second.re, (point.first + " re").c_str());
                                }
                                if (fabs(point.second.im) >= GMX_FLOAT_MIN)
                                {
                                    gridValuesChecker.checkReal(point.second.im, (point.first + " im").c_str());
                                }
                            }

                            if (computeEnergyAndVirial)
                            {
                                real       energy;
                                Matrix3x3  virial;
                                std::tie(energy, virial) = pmeGetReciprocalEnergyAndVirial(pmeSafe.get(), mode.first, method);
                                /* Energy */
                                TestReferenceChecker energyChecker(checker);
                                energyChecker.setDefaultTolerance(relativeToleranceAsUlp(10.0, 40));
                                energyChecker.checkReal(energy, "Energy");
                                /* Virial */
                                TestReferenceChecker virialChecker(checker.checkCompound("Matrix", "Virial"));
                                virialChecker.setDefaultTolerance(relativeToleranceAsUlp(1000, 8));
                                for (int i = 0; i < DIM; i++)
                                {
                                    for (int j = 0; j <= i; j++)
                                    {
                                        std::string valueId = formatString("Cell %d %d", i, j);
                                        virialChecker.checkReal(virial[i * DIM + j], valueId.c_str());
                                    }
                                }
                            }
                        }
                    }
                }

            }
        }
};

/*! \brief Test for PME solving */
TEST_P(PmeSolveTest, ReproducesOutputs)
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

//! A couple of valid inputs for grid sizes
static std::vector<IVec> const c_sampleGridSizes
{
    IVec {
        16, 12, 28
    },
    IVec {
        9, 7, 23
    }
};

//! Moved out from instantiations for readability
const auto c_inputBoxes     = ::testing::ValuesIn(c_sampleBoxes);
//! Moved out from instantiations for readability
const auto c_inputGridSizes = ::testing::ValuesIn(c_sampleGridSizes);

//! 2 sample complex grids - only non-zero values have to be listed
static std::vector<SparseComplexGridValuesInput> const c_sampleGrids
{
    SparseComplexGridValuesInput {{
                                      IVec {
                                          0, 0, 0
                                      }, t_complex {
                                          3.5f, 6.7f
                                      }
                                  }, {
                                      IVec {
                                          7, 0, 0
                                      }, t_complex {
                                          -2.5f, -0.7f
                                      }
                                  }, {
                                      IVec {
                                          3, 5, 7
                                      }, t_complex {
                                          -0.006f, 1e-8f
                                      }
                                  }, {
                                      IVec {
                                          3, 1, 2
                                      }, t_complex {
                                          0.6f, 7.9f
                                      }
                                  },  {
                                      IVec {
                                          6, 2, 4
                                      }, t_complex {
                                          30.1f, 2.45f
                                      }
                                  }, },
    SparseComplexGridValuesInput {{
                                      IVec {
                                          0, 4, 0
                                      }, t_complex {
                                          0.0f, 0.3f
                                      }
                                  }, {
                                      IVec {
                                          4, 2, 7
                                      }, t_complex {
                                          13.76f, -40.0f
                                      }
                                  }, {
                                      IVec {
                                          0, 6, 7
                                      }, t_complex {
                                          3.6f, 0.0f
                                      }
                                  }, {
                                      IVec {
                                          2, 5, 10
                                      }, t_complex {
                                          3.6f, 10.65f
                                      }
                                  }, }
};

//! Moved out from instantiations for readability
const auto c_inputGrids = ::testing::ValuesIn(c_sampleGrids);
//! Moved out from instantiations for readability
const auto c_inputEpsilon_r = ::testing::Values(1.2);
//! Moved out from instantiations for readability
const auto c_inputEwaldCoeff_q = ::testing::Values(2.0);
//! Moved out from instantiations for readability
const auto c_inputEwaldCoeff_lj = ::testing::Values(0.7);
//! Moved out from instantiations for readability
const auto c_inputMethods = ::testing::Values(PmeSolveAlgorithm::Coulomb, PmeSolveAlgorithm::LennardJones);

//! Instantiation of the PME solving test
INSTANTIATE_TEST_CASE_P(SaneInput, PmeSolveTest, ::testing::Combine(c_inputBoxes, c_inputGridSizes, c_inputGrids,
                                                                    c_inputEpsilon_r, c_inputEwaldCoeff_q, c_inputEwaldCoeff_lj, c_inputMethods));

//! A few more instances to check that different ewaldCoeff_q actually affects results of the Coulomb solver
INSTANTIATE_TEST_CASE_P(DifferentEwaldCoeffQ, PmeSolveTest, ::testing::Combine(c_inputBoxes, c_inputGridSizes, c_inputGrids,
                                                                               c_inputEpsilon_r, ::testing::Values(0.4), c_inputEwaldCoeff_lj,
                                                                                   ::testing::Values(PmeSolveAlgorithm::Coulomb)));

//! A few more instances to check that different ewaldCoeff_lj actually affects results of the Lennard-Jones solver.
//! The value has to be approximately larger than 1 / (box dimensions) to have a meaningful output grid.
//! Previous value of 0.3 caused one of the grid cells to be less or greater than GMX_FLOAT_MIN, depending on the architecture.
INSTANTIATE_TEST_CASE_P(DifferentEwaldCoeffLJ, PmeSolveTest, ::testing::Combine(c_inputBoxes, c_inputGridSizes, c_inputGrids,
                                                                                c_inputEpsilon_r, c_inputEwaldCoeff_q, ::testing::Values(2.35),
                                                                                    ::testing::Values(PmeSolveAlgorithm::LennardJones)));

//! A few more instances to check that different epsilon_r actually affects results of all solvers
INSTANTIATE_TEST_CASE_P(DifferentEpsilonR, PmeSolveTest, ::testing::Combine(c_inputBoxes, c_inputGridSizes, c_inputGrids,
                                                                            testing::Values(1.9), c_inputEwaldCoeff_q, c_inputEwaldCoeff_lj,
                                                                            c_inputMethods));

}
}
}
