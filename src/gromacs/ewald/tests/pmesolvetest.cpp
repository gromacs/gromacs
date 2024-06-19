/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * \brief
 * Implements PME solving tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include <cmath>
#include <cstdint>

#include <algorithm>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_gpu_internal.h"
#include "gromacs/ewald/pme_output.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"
#include "testutils/testinit.h"

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

/*! \brief Convenience typedef of the test input parameters
 *
 * Parameters:
 * - unit cell box
 * - complex grid dimensions
 * - complex grid values
 * - electrostatic constant epsilon_r
 * - Ewald splitting parameter ewaldcoeff_q
 * - Ewald splitting parameter ewaldcoeff_lj
 * - solver type
 * - grid ordering
 * - whether to compute energy and virial
 * - PME hardware context index
 *
 * Output: transformed local grid (Fourier space); optionally reciprocal energy and virial matrix.
 * TODO:
 * Implement and test Lorentz-Berthelot
 */
typedef std::tuple<std::string, IVec, std::string, double, double, double, PmeSolveAlgorithm, GridOrdering, bool, int> SolveInputParameters;

const char* enumValueToString(PmeSolveAlgorithm enumValue)
{
    static constexpr gmx::EnumerationArray<PmeSolveAlgorithm, const char*> s_strings = { "Coulomb",
                                                                                         "LJ" };
    return s_strings[enumValue];
}

const char* enumValueToString(GridOrdering enumValue)
{
    static constexpr gmx::EnumerationArray<GridOrdering, const char*> s_strings = {
        "YZX",
        "XYZ",
    };
    return s_strings[enumValue];
}

/*! \brief Help GoogleTest name our test cases
 *
 * This is intended to work like a custom test-naming function that
 * would be passed as the fourth argument to INSTANTIATE_TEST_SUITE_P,
 * except that we are not using that macro for these tests. Only the
 * components of SolveInputParameters that affect the reference data
 * values affect this name. Hardware context, grid ordering, and
 * whether this test targets energy&virial computation do not affect
 * this name. */
std::string nameOfTest(const testing::TestParamInfo<SolveInputParameters>& info)
{
    std::string testName = formatString(
            "box_%s_"
            "grid_%d_%d_%d_"
            "gridvalues_%s_"
            "eps_%g_"
            "ewaldq_%g_"
            "ewaldlj_%g_"
            "%s",
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

/*! \brief Help GoogleTest name our test cases
 *
 * This is intended to work like a custom test-naming function that
 * would be passed as the fourth argument to INSTANTIATE_TEST_SUITE_P,
 * except that we are not using that macro for these tests. All
 * components of SolveInputParameters affect this name. */
std::string fullNameOfTest(const testing::TestParamInfo<SolveInputParameters>& info,
                           const std::string&                                  testName)
{
    // Note that makeRefDataFileName() relies on finding "WorksOn" in
    // the name of the test case, so it can remove the information
    // about the hardware context and all following text from the name
    // of the file used for refdata.
    const int hardwareContextIndex = std::get<9>(info.param);
    return formatString(
            "WorksOn%s_%s_order_%s_"
            "%s",
            makeHardwareContextName(hardwareContextIndex).c_str(),
            testName.c_str(),
            enumValueToString(std::get<7>(info.param)),
            std::get<8>(info.param) ? "energy" : "");
}

//! Test fixture
class SolveTest : public ::testing::TestWithParam<SolveInputParameters>
{
public:
    SolveTest() = default;
};

/*! \brief Test case whose body checks that spline and spread work
 *
 * Normally the declaration of this class would be produced by a call
 * to a macro like TEST_P(SolveTest, WorksWith). That macro places the
 * body of the test case in the TestBody() method, which here is done
 * explicitly.
 *
 * Note that it is important to use parameters_ to access the values
 * that describe the particular test case, rather than the usual
 * GoogleTest function GetParam(), because the latter no longer
 * works. */
class SolveTestBody : public SolveTest
{
public:
    //! Constructor
    explicit SolveTestBody(const SolveInputParameters& parameters) : parameters_(parameters) {}

    //! The test parameters with which the test case was instantiated
    SolveInputParameters parameters_;
    //! The test
    void TestBody() override
    {
        /* Getting the input */
        IVec              gridSize;
        double            epsilon_r;
        double            ewaldCoeff_q;
        double            ewaldCoeff_lj;
        PmeSolveAlgorithm method;
        std::string       boxName, gridValuesName;
        GridOrdering      gridOrdering;
        bool              computeEnergyAndVirial;
        int               contextIndex;
        std::tie(boxName,
                 gridSize,
                 gridValuesName,
                 epsilon_r,
                 ewaldCoeff_q,
                 ewaldCoeff_lj,
                 method,
                 gridOrdering,
                 computeEnergyAndVirial,
                 contextIndex)                                = parameters_;
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

        const auto                    pmeTestHardwareContextsAll = getPmeTestHardwareContexts();
        const PmeTestHardwareContext& pmeTestHardwareContext = pmeTestHardwareContextsAll[contextIndex];
        const CodePath                codePath               = pmeTestHardwareContext.codePath();
        MessageStringCollector        messages = getSkipMessagesIfNecessary(inputRec, codePath);
        messages.appendIf(!pmeTestHardwareContext.gpuId().has_value() && gridOrdering == GridOrdering::XYZ,
                          "CPU PME solve does not implement XYZ grid ordering");
        if (!messages.isEmpty())
        {
            GTEST_SKIP() << messages.toString();
        }

        pmeTestHardwareContext.activate();
        SCOPED_TRACE("Testing on " + pmeTestHardwareContext.description());

        /* Describing the test*/
        SCOPED_TRACE(formatString(
                "Testing solving (%s, %s, %s energy/virial) on %s for PME grid "
                "size %d %d %d, Ewald coefficients %g %g",
                (method == PmeSolveAlgorithm::LennardJones) ? "Lennard-Jones" : "Coulomb",
                enumValueToString(gridOrdering),
                computeEnergyAndVirial ? "with" : "without",
                pmeTestHardwareContext.description().c_str(),
                gridSize[XX],
                gridSize[YY],
                gridSize[ZZ],
                ewaldCoeff_q,
                ewaldCoeff_lj));

        /* Running the test */
        PmeSafePointer pmeSafe = pmeInitWrapper(&inputRec,
                                                codePath,
                                                pmeTestHardwareContext.deviceContext(),
                                                pmeTestHardwareContext.deviceStream(),
                                                pmeTestHardwareContext.pmeGpuProgram(),
                                                box,
                                                ewaldCoeff_q,
                                                ewaldCoeff_lj);
        pmeSetComplexGrid(pmeSafe.get(), codePath, gridOrdering, nonZeroGridValues);
        const real cellVolume = box[0] * box[4] * box[8];
        // FIXME - this is box[XX][XX] * box[YY][YY] * box[ZZ][ZZ], should be stored in the PME structure
        pmePerformSolve(pmeSafe.get(), codePath, method, cellVolume, gridOrdering, computeEnergyAndVirial);
        pmeFinalizeTest(pmeSafe.get(), codePath);

        /* Check the outputs */
        TestReferenceData    refData(makeRefDataFileName());
        TestReferenceChecker checker(refData.rootChecker());

        SparseComplexGridValuesOutput nonZeroGridValuesOutput =
                pmeGetComplexGrid(pmeSafe.get(), codePath, gridOrdering);
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
            if (std::fabs(point.second.re) >= GMX_FLOAT_MIN)
            {
                gridValuesChecker.checkReal(point.second.re, (point.first + " re").c_str());
            }
            if (std::fabs(point.second.im) >= GMX_FLOAT_MIN)
            {
                gridValuesChecker.checkReal(point.second.im, (point.first + " im").c_str());
            }
        }

        TestReferenceChecker energyChecker(checker);
        TestReferenceChecker virialChecker(checker.checkCompound("Matrix", "Virial"));
        if (computeEnergyAndVirial)
        {
            // Extract the energy and virial
            const auto  output = pmeGetReciprocalEnergyAndVirial(pmeSafe.get(), codePath, method);
            const auto& energy = (method == PmeSolveAlgorithm::Coulomb) ? output.coulombEnergy_
                                                                        : output.lennardJonesEnergy_;
            const auto& virial = (method == PmeSolveAlgorithm::Coulomb) ? output.coulombVirial_
                                                                        : output.lennardJonesVirial_;

            // These quantities are computed based on the grid values, so must have
            // checking relative tolerances at least as large. Virial needs more flops
            // than energy, so needs a larger tolerance.

            /* Energy */
            double energyMagnitude = 10.0;
            // TODO This factor is arbitrary, do a proper error-propagation analysis
            uint64_t energyUlpToleranceFactor = gridUlpToleranceFactor * 2;
            auto     energyTolerance          = relativeToleranceAsPrecisionDependentUlp(
                    energyMagnitude,
                    energyUlpToleranceFactor * c_splineModuliSinglePrecisionUlps,
                    energyUlpToleranceFactor * splineModuliDoublePrecisionUlps);
            energyChecker.setDefaultTolerance(energyTolerance);
            energyChecker.checkReal(energy, "Energy");

            /* Virial */
            double virialMagnitude = 1000.0;
            // TODO This factor is arbitrary, do a proper error-propagation analysis
            uint64_t virialUlpToleranceFactor = energyUlpToleranceFactor * 2;
            auto     virialTolerance          = relativeToleranceAsPrecisionDependentUlp(
                    virialMagnitude,
                    virialUlpToleranceFactor * c_splineModuliSinglePrecisionUlps,
                    virialUlpToleranceFactor * splineModuliDoublePrecisionUlps);
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
        else
        {
            energyChecker.disableUnusedEntriesCheck();
            virialChecker.disableUnusedEntriesCheck();
        }
    }
};

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
//! Moved out from instantiations for readability
const auto c_gridOrderings = ::testing::Values(GridOrdering::XYZ, GridOrdering::YZX);

} // namespace

void registerDynamicalPmeSolveTests(const Range<int> hardwareContextIndexRange)
{
    {
        // Form the Cartesian product of all test values we might check
        auto testCombinations = ::testing::Combine(
                c_inputBoxNames,
                ::testing::ValuesIn(c_inputGridSizes),
                c_inputGridNames,
                c_inputEpsilon_r,
                c_inputEwaldCoeff_q,
                c_inputEwaldCoeff_lj,
                c_inputMethods,
                c_gridOrderings,
                ::testing::Bool(),
                ::testing::Range(*hardwareContextIndexRange.begin(), *hardwareContextIndexRange.end()));
        gmx::test::registerTests<SolveTest, SolveTestBody, decltype(testCombinations)>(
                "Pme_SolveTest", nameOfTest, fullNameOfTest, testCombinations);
    }
    {
        // A few more instances to check that different ewaldCoeff_q actually affects results of the Coulomb solver
        auto testCombinations = ::testing::Combine(
                c_inputBoxNames,
                ::testing::ValuesIn(c_inputGridSizes),
                c_inputGridNames,
                c_inputEpsilon_r,
                ::testing::Values(0.4),
                c_inputEwaldCoeff_lj,
                ::testing::Values(PmeSolveAlgorithm::Coulomb),
                c_gridOrderings,
                ::testing::Bool(),
                ::testing::Range(*hardwareContextIndexRange.begin(), *hardwareContextIndexRange.end()));
        gmx::test::registerTests<SolveTest, SolveTestBody, decltype(testCombinations)>(
                "PmeDiffEwaldQ_SolveTest", nameOfTest, fullNameOfTest, testCombinations);
    }
    {
        // A few more instances to check that different ewaldCoeff_lj
        // actually affects results of the Lennard-Jones solver.  The
        // value has to be approximately larger than 1 / (box dimensions)
        // to have a meaningful output grid.  Previous value of 0.3 caused
        // one of the grid cells to be less or greater than GMX_FLOAT_MIN,
        // depending on the architecture.
        auto testCombinations = ::testing::Combine(
                c_inputBoxNames,
                ::testing::ValuesIn(c_inputGridSizes),
                c_inputGridNames,
                c_inputEpsilon_r,
                c_inputEwaldCoeff_q,
                ::testing::Values(2.35),
                ::testing::Values(PmeSolveAlgorithm::LennardJones),
                c_gridOrderings,
                ::testing::Bool(),
                ::testing::Range(*hardwareContextIndexRange.begin(), *hardwareContextIndexRange.end()));
        gmx::test::registerTests<SolveTest, SolveTestBody, decltype(testCombinations)>(
                "PmeDiffEwaldLJ_SolveTest", nameOfTest, fullNameOfTest, testCombinations);
    }
    {
        // A few more instances to check that different epsilon_r actually affects results of all solvers
        auto testCombinations = ::testing::Combine(
                c_inputBoxNames,
                ::testing::ValuesIn(c_inputGridSizes),
                c_inputGridNames,
                testing::Values(1.9),
                c_inputEwaldCoeff_q,
                c_inputEwaldCoeff_lj,
                c_inputMethods,
                c_gridOrderings,
                ::testing::Bool(),
                ::testing::Range(*hardwareContextIndexRange.begin(), *hardwareContextIndexRange.end()));
        gmx::test::registerTests<SolveTest, SolveTestBody, decltype(testCombinations)>(
                "PmeDiffEps_SolveTest", nameOfTest, fullNameOfTest, testCombinations);
    }
}

} // namespace test
} // namespace gmx
