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
/*! \brief Convenience typedef of the test input parameters - unit cell box, complex grid dimensions, complex grid values,
 * electrostatic constant epsilon_r, Ewald splitting parameters ewaldcoeff_q and ewaldcoeff_lj, solver type
 * Output: transformed local grid (Fourier space); optionally reciprocal energy and virial matrix.
 * TODO:
 * Test bEnerVir off as well? Test Lorentz-Berthelot? Test XYZ / YZX for GPU case?
 */

typedef std::tuple<Matrix3x3, IVec, SparseComplexGridValuesInput, real, real, real, PmeSolveAlgorithm> SolveInputParameters;

//! Test fixture
class PmeSolveTest : public ::testing::TestWithParam<SolveInputParameters>
{
    private:
        //! Environment for getting the t_inputrec structure easily
        MDModules mdModules_;

    public:
        //! Default constructor
        PmeSolveTest() = default;

        //! The test
        void runTest()
        {
            /* Getting the input */
            Matrix3x3                    box;
            IVec                         gridSize;
            SparseComplexGridValuesInput nonZeroGridValues;
            real                         epsilon_r;
            real                         ewaldCoeff_q;
            real                         ewaldCoeff_lj;
            PmeSolveAlgorithm            method;
            std::tie(box, gridSize, nonZeroGridValues, epsilon_r, ewaldCoeff_q, ewaldCoeff_lj, method) = GetParam();

            /* Storing the input where it's needed, running the test */
            t_inputrec *inputRec  = mdModules_.inputrec();
            inputRec->nkx         = gridSize[XX];
            inputRec->nky         = gridSize[YY];
            inputRec->nkz         = gridSize[ZZ];
            inputRec->pme_order   = 4;
            inputRec->coulombtype = eelPME;
            inputRec->epsilon_r   = epsilon_r;
            switch (method)
            {
                case PmeSolveAlgorithm::Normal:
                    break;

                case PmeSolveAlgorithm::LennardJones:
                    inputRec->vdwtype = evdwPME;
                    break;

                default:
                    GMX_THROW(InternalError("Unknown PME solver"));
            }

            TestReferenceData                     refData;
            const std::map<CodePath, std::string> modesToTest = {{CodePath::CPU, "CPU"}};
            for (const auto &mode : modesToTest)
            {
                /* Describing the test*/
                SCOPED_TRACE(formatString("Testing solving (%s) with %s for PME grid size %d %d %d, Ewald coefficients %g %g",
                                          (method == PmeSolveAlgorithm::LennardJones) ? "Lennard-Jones" : "Normal",
                                          mode.second.c_str(),
                                          gridSize[XX], gridSize[YY], gridSize[ZZ],
                                          ewaldCoeff_q, ewaldCoeff_lj
                                          ));

                /* Running the test */
                PmeSafePointer pmeSafe = pmeInitEmpty(inputRec, mode.first, box, ewaldCoeff_q, ewaldCoeff_lj);
                pmeSetComplexGrid(pmeSafe.get(), mode.first, nonZeroGridValues);
                const real     cellVolume = box[0] * box[4] * box[8];
                //FIXME - this is box[XX][XX] * box[YY][YY] * box[ZZ][ZZ], should be stored in the PME structure
                pmePerformSolve(pmeSafe.get(), mode.first, method, cellVolume);

                /* Check the outputs */
                TestReferenceChecker checker(refData.rootChecker());
                const auto           ulpTolerance = 40;
                checker.setDefaultTolerance(relativeToleranceAsUlp(1.0, ulpTolerance));

                SparseComplexGridValuesOutput nonZeroGridValuesOutput = pmeGetComplexGrid(pmeSafe.get(), mode.first);
                /* Transformed grid */
                TestReferenceChecker          gridValuesChecker(checker.checkCompound("NonZeroGridValues", "ComplexSpaceGrid"));
                const auto                    ulpToleranceGrid = 4;
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

                real       energy;
                Matrix3x3  virial;
                std::tie(energy, virial) = pmeGetReciprocalEnergyAndVirial(pmeSafe.get(), mode.first, method);
                /* Energy */
                checker.checkReal(energy, "Energy");
                /* Virial */
                TestReferenceChecker virialChecker(checker.checkCompound("Matrix", "Virial"));
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
};

/*! \brief Test for PME solving */
TEST_P(PmeSolveTest, ReproducesOutputs)
{
    EXPECT_NO_THROW(runTest());
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
                   0.0f, 3.1f, 0.0f,
                   6.0f, 4.0f, 3.2f
               }},
};

//! A couple of valid inputs for grid sizes
static std::vector<IVec> const sampleGridSizes
{
    IVec {
        16, 12, 28
    },
    IVec {
        9, 7, 23
    }
};

//! Moved out from instantiations for readability
auto inputBoxes     = ::testing::ValuesIn(sampleBoxes);
//! Moved out from instantiations for readability
auto inputGridSizes = ::testing::ValuesIn(sampleGridSizes);

//! 2 sample complex grids - only non-zero values have to be listed
static std::vector<SparseComplexGridValuesInput> const sampleGrids
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
auto inputGrids = ::testing::ValuesIn(sampleGrids);
//! Moved out from instantiations for readability
auto inputEpsilon_r = ::testing::Values(0.6, 1.0);
//! Moved out from instantiations for readability
auto inputEwaldCoeff_q = ::testing::Values(0.4, 2.0);
//! Moved out from instantiations for readability
auto inputEwaldCoeff_lj = ::testing::Values(0.7, 1.3);
//! Moved out from instantiations for readability
auto inputMethods = ::testing::Values(PmeSolveAlgorithm::Normal, PmeSolveAlgorithm::LennardJones);

//! Instantiation of the PME solving test
INSTANTIATE_TEST_CASE_P(SaneInput, PmeSolveTest, ::testing::Combine(inputBoxes, inputGridSizes, inputGrids,
                                                                    inputEpsilon_r, inputEwaldCoeff_q, inputEwaldCoeff_lj, inputMethods));
}
}
}
