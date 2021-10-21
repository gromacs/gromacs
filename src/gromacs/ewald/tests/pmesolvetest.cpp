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
 * Implements PME solving tests.
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
std::vector<IVec> const c_inputGridSizes{ IVec{ 16, 12, 28 }, IVec{ 9, 7, 23 } };

//! Two input complex grids - only non-zero values have to be listed
const std::map<std::string, SparseComplexGridValuesInput> c_inputGridValues = {
    { "first",
      SparseComplexGridValuesInput{
              { IVec{ 0, 0, 0 }, t_complex{ 3.5F, 6.7F } },
              { IVec{ 7, 0, 0 }, t_complex{ -2.5F, -0.7F } },
              { IVec{ 3, 5, 7 }, t_complex{ -0.006F, 1e-8F } },
              { IVec{ 3, 1, 2 }, t_complex{ 0.6F, 7.9F } },
              { IVec{ 6, 2, 4 }, t_complex{ 30.1F, 2.45F } },
      } },
    { "second",
      SparseComplexGridValuesInput{
              { IVec{ 0, 4, 0 }, t_complex{ 0.0F, 0.3F } },
              { IVec{ 4, 2, 7 }, t_complex{ 13.76F, -40.0F } },
              { IVec{ 0, 6, 7 }, t_complex{ 3.6F, 0.0F } },
              { IVec{ 2, 5, 10 }, t_complex{ 3.6F, 10.65F } },
      } },
};

/*! \brief Convenience typedef of the test input parameters - unit cell box, complex grid dimensions, complex grid values,
 * electrostatic constant epsilon_r, Ewald splitting parameters ewaldcoeff_q and ewaldcoeff_lj, solver type
 * Output: transformed local grid (Fourier space); optionally reciprocal energy and virial matrix.
 * TODO:
 * Implement and test Lorentz-Berthelot
 */
typedef std::tuple<std::string, IVec, std::string, double, double, double, PmeSolveAlgorithm> SolveInputParameters;

const char* enumValueToString(PmeSolveAlgorithm enumValue)
{
    static constexpr gmx::EnumerationArray<PmeSolveAlgorithm, const char*> s_pmeSolveAlgorithmNames = {
        "Coulomb", "LJ"
    };
    return s_pmeSolveAlgorithmNames[enumValue];
}

//! Help GoogleTest name our test cases
std::string nameOfTest(const testing::TestParamInfo<SolveInputParameters>& info)
{
    std::string testName = formatString(
            "box_%s_"
            "grid_%d_%d_%d_"
            "gridvalues_%s_"
            "eps_%g_"
            "ewaldq_%g_"
            "ewaldlj_%g_"
            "method_%s",
            std::get<0>(info.param).c_str(),
            std::get<1>(info.param)[XX],
            std::get<1>(info.param)[YY],
            std::get<1>(info.param)[ZZ],
            std::get<2>(info.param).c_str(),
            std::get<3>(info.param),
            std::get<4>(info.param),
            std::get<5>(info.param),
            enumValueToString(std::get<6>(info.param)));

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

//! Test fixture
class SolveTest : public ::testing::TestWithParam<SolveInputParameters>
{
public:
    SolveTest() = default;

    //! Sets the programs once
    static void SetUpTestSuite()
    {
        s_pmeTestHardwareContexts    = createPmeTestHardwareContextList();
        g_allowPmeWithSyclForTesting = true; // We support Solve with SYCL
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
        IVec              gridSize;
        double            epsilon_r;
        double            ewaldCoeff_q;
        double            ewaldCoeff_lj;
        PmeSolveAlgorithm method;
        std::string       boxName, gridValuesName;
        std::tie(boxName, gridSize, gridValuesName, epsilon_r, ewaldCoeff_q, ewaldCoeff_lj, method) =
                GetParam();
        Matrix3x3                           box               = c_inputBoxes.at(boxName);
        const SparseComplexGridValuesInput& nonZeroGridValues = c_inputGridValues.at(gridValuesName);

        /* Storing the input where it's needed, running the test */
        t_inputrec inputRec;
        inputRec.nkx         = gridSize[XX];
        inputRec.nky         = gridSize[YY];
        inputRec.nkz         = gridSize[ZZ];
        inputRec.pme_order   = 4;
        inputRec.coulombtype = CoulombInteractionType::Pme;
        inputRec.epsilon_r   = epsilon_r;
        switch (method)
        {
            case PmeSolveAlgorithm::Coulomb: break;

            case PmeSolveAlgorithm::LennardJones: inputRec.vdwtype = VanDerWaalsType::Pme; break;

            default: GMX_THROW(InternalError("Unknown PME solver"));
        }

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
                EXPECT_THROW_GMX(
                        pmeInitWrapper(&inputRec, codePath, nullptr, nullptr, nullptr, box, ewaldCoeff_q, ewaldCoeff_lj),
                        NotImplementedError);
                continue;
            }

            std::map<GridOrdering, std::string> gridOrderingsToTest = { { GridOrdering::YZX,
                                                                          "YZX" } };
            if (codePath == CodePath::GPU)
            {
                gridOrderingsToTest[GridOrdering::XYZ] = "XYZ";
            }
            for (const auto& gridOrdering : gridOrderingsToTest)
            {
                for (bool computeEnergyAndVirial : { false, true })
                {
                    /* Describing the test*/
                    SCOPED_TRACE(formatString(
                            "Testing solving (%s, %s, %s energy/virial) on %s for PME grid "
                            "size %d %d %d, Ewald coefficients %g %g",
                            (method == PmeSolveAlgorithm::LennardJones) ? "Lennard-Jones" : "Coulomb",
                            gridOrdering.second.c_str(),
                            computeEnergyAndVirial ? "with" : "without",
                            pmeTestHardwareContext->description().c_str(),
                            gridSize[XX],
                            gridSize[YY],
                            gridSize[ZZ],
                            ewaldCoeff_q,
                            ewaldCoeff_lj));

                    /* Running the test */
                    PmeSafePointer pmeSafe = pmeInitWrapper(&inputRec,
                                                            codePath,
                                                            pmeTestHardwareContext->deviceContext(),
                                                            pmeTestHardwareContext->deviceStream(),
                                                            pmeTestHardwareContext->pmeGpuProgram(),
                                                            box,
                                                            ewaldCoeff_q,
                                                            ewaldCoeff_lj);
                    pmeSetComplexGrid(pmeSafe.get(), codePath, gridOrdering.first, nonZeroGridValues);
                    const real cellVolume = box[0] * box[4] * box[8];
                    // FIXME - this is box[XX][XX] * box[YY][YY] * box[ZZ][ZZ], should be stored in the PME structure
                    pmePerformSolve(pmeSafe.get(), codePath, method, cellVolume, gridOrdering.first, computeEnergyAndVirial);
                    pmeFinalizeTest(pmeSafe.get(), codePath);

                    /* Check the outputs */
                    TestReferenceChecker checker(refData.rootChecker());

                    SparseComplexGridValuesOutput nonZeroGridValuesOutput =
                            pmeGetComplexGrid(pmeSafe.get(), codePath, gridOrdering.first);
                    /* Transformed grid */
                    TestReferenceChecker gridValuesChecker(
                            checker.checkCompound("NonZeroGridValues", "ComplexSpaceGrid"));

                    real gridValuesMagnitude = 1.0;
                    for (const auto& point : nonZeroGridValuesOutput)
                    {
                        gridValuesMagnitude = std::max(std::fabs(point.second.re), gridValuesMagnitude);
                        gridValuesMagnitude = std::max(std::fabs(point.second.im), gridValuesMagnitude);
                    }
                    // Spline moduli participate 3 times in the computation; 2 is an additional factor for SIMD exp() precision
                    uint64_t gridUlpToleranceFactor = DIM * 2;
                    if (method == PmeSolveAlgorithm::LennardJones)
                    {
                        // Lennard Jones is more complex and also uses erfc(), relax more
                        gridUlpToleranceFactor *= 2;
                    }
                    const uint64_t splineModuliDoublePrecisionUlps =
                            getSplineModuliDoublePrecisionUlps(inputRec.pme_order + 1);
                    auto gridTolerance = relativeToleranceAsPrecisionDependentUlp(
                            gridValuesMagnitude,
                            gridUlpToleranceFactor * c_splineModuliSinglePrecisionUlps,
                            gridUlpToleranceFactor * splineModuliDoublePrecisionUlps);
                    gridValuesChecker.setDefaultTolerance(gridTolerance);

                    for (const auto& point : nonZeroGridValuesOutput)
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
                        // Extract the energy and virial
                        const auto output =
                                pmeGetReciprocalEnergyAndVirial(pmeSafe.get(), codePath, method);
                        const auto& energy = (method == PmeSolveAlgorithm::Coulomb)
                                                     ? output.coulombEnergy_
                                                     : output.lennardJonesEnergy_;
                        const auto& virial = (method == PmeSolveAlgorithm::Coulomb)
                                                     ? output.coulombVirial_
                                                     : output.lennardJonesVirial_;

                        // These quantities are computed based on the grid values, so must have
                        // checking relative tolerances at least as large. Virial needs more flops
                        // than energy, so needs a larger tolerance.

                        /* Energy */
                        double energyMagnitude = 10.0;
                        // TODO This factor is arbitrary, do a proper error-propagation analysis
                        uint64_t energyUlpToleranceFactor = gridUlpToleranceFactor * 2;
                        auto     energyTolerance = relativeToleranceAsPrecisionDependentUlp(
                                energyMagnitude,
                                energyUlpToleranceFactor * c_splineModuliSinglePrecisionUlps,
                                energyUlpToleranceFactor * splineModuliDoublePrecisionUlps);
                        TestReferenceChecker energyChecker(checker);
                        energyChecker.setDefaultTolerance(energyTolerance);
                        energyChecker.checkReal(energy, "Energy");

                        /* Virial */
                        double virialMagnitude = 1000.0;
                        // TODO This factor is arbitrary, do a proper error-propagation analysis
                        uint64_t virialUlpToleranceFactor = energyUlpToleranceFactor * 2;
                        auto     virialTolerance = relativeToleranceAsPrecisionDependentUlp(
                                virialMagnitude,
                                virialUlpToleranceFactor * c_splineModuliSinglePrecisionUlps,
                                virialUlpToleranceFactor * splineModuliDoublePrecisionUlps);
                        TestReferenceChecker virialChecker(
                                checker.checkCompound("Matrix", "Virial"));
                        virialChecker.setDefaultTolerance(virialTolerance);
                        for (int i = 0; i < DIM; i++)
                        {
                            for (int j = 0; j <= i; j++)
                            {
                                std::string valueId = formatString("Cell %d %d", i, j);
                                virialChecker.checkReal(virial[i][j], valueId.c_str());
                            }
                        }
                    }
                }
            }
        }
    }

    static std::vector<std::unique_ptr<PmeTestHardwareContext>> s_pmeTestHardwareContexts;
};

std::vector<std::unique_ptr<PmeTestHardwareContext>> SolveTest::s_pmeTestHardwareContexts;

/*! \brief Test for PME solving */
TEST_P(SolveTest, WorksWith)
{
    EXPECT_NO_THROW_GMX(runTest());
}

//! Moved out from instantiations for readability
const auto c_inputBoxNames = ::testing::Values("rect", "tric");
//! Moved out from instantiations for readability
const auto c_inputGridNames = ::testing::Values("first", "second");
//! Moved out from instantiations for readability
const auto c_inputEpsilon_r = ::testing::Values(1.2);
//! Moved out from instantiations for readability
const auto c_inputEwaldCoeff_q = ::testing::Values(2.0);
//! Moved out from instantiations for readability
const auto c_inputEwaldCoeff_lj = ::testing::Values(0.7);
//! Moved out from instantiations for readability
const auto c_inputMethods = ::testing::Values(PmeSolveAlgorithm::Coulomb, PmeSolveAlgorithm::LennardJones);

//! Instantiation of the PME solving test
INSTANTIATE_TEST_SUITE_P(Pme,
                         SolveTest,
                         ::testing::Combine(c_inputBoxNames,
                                            ::testing::ValuesIn(c_inputGridSizes),
                                            c_inputGridNames,
                                            c_inputEpsilon_r,
                                            c_inputEwaldCoeff_q,
                                            c_inputEwaldCoeff_lj,
                                            c_inputMethods),
                         nameOfTest);

//! A few more instances to check that different ewaldCoeff_q actually affects results of the Coulomb solver
INSTANTIATE_TEST_SUITE_P(PmeDiffEwaldQ,
                         SolveTest,
                         ::testing::Combine(c_inputBoxNames,
                                            ::testing::ValuesIn(c_inputGridSizes),
                                            c_inputGridNames,
                                            c_inputEpsilon_r,
                                            ::testing::Values(0.4),
                                            c_inputEwaldCoeff_lj,
                                            ::testing::Values(PmeSolveAlgorithm::Coulomb)),
                         nameOfTest);

//! A few more instances to check that different ewaldCoeff_lj actually affects results of the Lennard-Jones solver.
//! The value has to be approximately larger than 1 / (box dimensions) to have a meaningful output grid.
//! Previous value of 0.3 caused one of the grid cells to be less or greater than GMX_FLOAT_MIN, depending on the architecture.
INSTANTIATE_TEST_SUITE_P(PmeDiffEwaldLJ,
                         SolveTest,
                         ::testing::Combine(c_inputBoxNames,
                                            ::testing::ValuesIn(c_inputGridSizes),
                                            c_inputGridNames,
                                            c_inputEpsilon_r,
                                            c_inputEwaldCoeff_q,
                                            ::testing::Values(2.35),
                                            ::testing::Values(PmeSolveAlgorithm::LennardJones)),
                         nameOfTest);

//! A few more instances to check that different epsilon_r actually affects results of all solvers
INSTANTIATE_TEST_SUITE_P(PmeDiffEps,
                         SolveTest,
                         ::testing::Combine(c_inputBoxNames,
                                            ::testing::ValuesIn(c_inputGridSizes),
                                            c_inputGridNames,
                                            testing::Values(1.9),
                                            c_inputEwaldCoeff_q,
                                            c_inputEwaldCoeff_lj,
                                            c_inputMethods),
                         nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
