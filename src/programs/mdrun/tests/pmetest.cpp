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
 * This implements basic PME sanity tests for end-to-end mdrun simulations.
 * It runs the input system with PME for several steps (on CPU and GPU, if available),
 * and checks the reciprocal and conserved energies.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdlib>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/mpiinfo.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! Enum describing the flavors of PME tests that are run
enum class PmeTestFlavor : int
{
    Basic,
    WithWalls,
    Count
};

//! Helper to print a string to describe the PME test flavor
const char* enumValueToString(const PmeTestFlavor enumValue)
{
    static constexpr gmx::EnumerationArray<PmeTestFlavor, const char*> s_names = {
        "basic",
        "with walls",
    };
    return s_names[enumValue];
}

// Paramters for parametrized test fixture: the flavor of PME test to
// run, and options for an mdrun command line.
using PmeTestParameters = std::tuple<PmeTestFlavor, std::string>;

/*! \brief Help GoogleTest name our tests
 *
 * If changes are needed here, consider making matching changes in
 * makeRefDataFileName(). */
std::string nameOfTest(const testing::TestParamInfo<PmeTestParameters>& info)
{
    std::string testName = formatString(
            "%s_mdrun_%s", enumValueToString(std::get<0>(info.param)), std::get<1>(info.param).c_str());

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "");
    testName = replaceAll(testName, " ", "_");
    return testName;
}

/*! \brief Construct a refdata filename for this test
 *
 * We want the same reference data to apply to every mdrun command
 * line that we test. That means we need to store it in a file whose
 * name relates to the name of the test excluding the part related to
 * the mdrun command line. By default, the reference data filename is
 * set via a call to gmx::TestFileManager::getTestSpecificFileName()
 * that queries GoogleTest and gets a string that includes the return
 * value for nameOfTest(). This code works similarly, but removes the
 * aforementioned part. This logic must match the implementation of
 * nameOfTest() so that it works as intended. */
std::string makeRefDataFileName()
{
    // Get the info about the test
    const ::testing::TestInfo* testInfo = ::testing::UnitTest::GetInstance()->current_test_info();

    // Get the test name and edit it to remove the mdrun command-line
    // part.
    std::string testName(testInfo->name());
    auto        separatorPos = testName.find("_mdrun");
    testName                 = testName.substr(0, separatorPos);

    // Build the complete filename like getTestSpecificFilename() does
    // it.
    std::string testSuiteName(testInfo->test_suite_name());
    std::string refDataFileName = testSuiteName + "_" + testName + ".xml";
    std::replace(refDataFileName.begin(), refDataFileName.end(), '/', '_');

    return refDataFileName;
}

/*! \brief Test fixture for end-to-end execution of PME */
class PmeTest : public MdrunTestFixture, public ::testing::WithParamInterface<PmeTestParameters>
{
public:
    //! Names of tpr files built by grompp in SetUpTestSuite to run in tests
    inline static EnumerationArray<PmeTestFlavor, std::string> s_tprFileNames;
    //! Mutex to protect creation of the TestFileManager with thread-MPI
    inline static std::mutex s_mutexForTestFileManager;
    //! Manager needed in SetUpTestSuite to handle making tpr files
    inline static std::unique_ptr<TestFileManager> s_testFileManager;
    //! Name of the system from the simulation database to use in SetupTestSuite
    inline static const std::string s_inputFile = "spc-and-methanol";
    //! Runs grompp to prepare the tpr files that are reused in the tests
    static void SetUpTestSuite();
    //! Cleans up the tpr files
    static void TearDownTestSuite();
    /*! \brief If the mdrun command line can't run on this build or
     * hardware, we should mark it as skipped and describe why. */
    static MessageStringCollector getSkipMessagesIfNecessary(const CommandLine& commandLine);
    //! Check the energies against the reference data.
    void checkEnergies(bool usePmeTuning) const;
};

// static
void PmeTest::SetUpTestSuite()
{
    MdrunTestFixture::SetUpTestSuite();
    // Ensure we only make one TestFileManager per process, which
    // ensures there is no race with thread-MPI. Whichever thread gets
    // the lock first initializes s_testFileManager, and by the time
    // the rest get the lock it is already initialized.  Without
    // thread-MPI, there's only one thread, so the initialization is
    // trivial.
    {
        std::lock_guard<std::mutex> testFileManagerLock(s_mutexForTestFileManager);
        if (!s_testFileManager)
        {
            s_testFileManager = std::make_unique<TestFileManager>();
        }
    }

    const static std::unordered_map<PmeTestFlavor, std::string> sc_pmeTestFlavorExtraMdpLines = {
        { PmeTestFlavor::Basic,
          { "nsteps = 20\n"
            "pbc    = xyz\n" } },
        { PmeTestFlavor::WithWalls,
          { "nsteps          = 0\n"
            "pbc             = xy\n"
            "nwall           = 2\n"
            "ewald-geometry  = 3dc\n"
            "wall_atomtype   = CMet H\n"
            "wall_density    = 9 9.0\n"
            "wall-ewald-zfac = 5\n" } },
    };

    // Make the tpr files for the different flavors
    for (PmeTestFlavor pmeTestFlavor : EnumerationWrapper<PmeTestFlavor>{})
    {
        SimulationRunner runner(s_testFileManager.get());
        runner.useTopGroAndNdxFromDatabase(s_inputFile);

        std::string mdpInputFileContents(
                "coulombtype     = PME\n"
                "nstcalcenergy   = 1\n"
                "nstenergy       = 1\n"
                "pme-order       = 4\n");
        mdpInputFileContents += sc_pmeTestFlavorExtraMdpLines.at(pmeTestFlavor);
        runner.useStringAsMdpFile(mdpInputFileContents);

        std::string tprFileNameSuffix = formatString("%s.tpr", enumValueToString(pmeTestFlavor));
        std::replace(tprFileNameSuffix.begin(), tprFileNameSuffix.end(), ' ', '_');
        runner.tprFileName_ = s_testFileManager->getTemporaryFilePath(tprFileNameSuffix).string();
        // Note that only one rank actually generates a tpr file
        runner.callGrompp();
        s_tprFileNames[pmeTestFlavor] = runner.tprFileName_;
    }
}

// static
void PmeTest::TearDownTestSuite()
{
    MdrunTestFixture::TearDownTestSuite();

    // Ensure we only clean up the TestFileManager once per process, which
    // ensures there is no race with thread-MPI.
    {
        std::lock_guard<std::mutex> testFileManagerLock(s_mutexForTestFileManager);
        s_testFileManager.reset(nullptr);
    }
}

MessageStringCollector PmeTest::getSkipMessagesIfNecessary(const CommandLine& commandLine)
{
    // Note that we can't call GTEST_SKIP() from within this method,
    // because it only returns from the current function. So we
    // collect all the reasons why the test cannot run, return them
    // and skip in a higher stack frame.

    MessageStringCollector messages;
    messages.startContext("Test is being skipped because:");

    const int numRanks = getNumberOfTestMpiRanks();

    // -npme was already required
    std::string npmeOptionArgument(commandLine.argumentOf("-npme").value());
    const bool  commandLineTargetsPmeOnlyRanks = (std::stoi(npmeOptionArgument) > 0);
    messages.appendIf(commandLineTargetsPmeOnlyRanks && numRanks == 1,
                      "it targets using PME rank(s) but the simulation is using only one rank");

    // -pme was already required
    std::string pmeOptionArgument(commandLine.argumentOf("-pme").value());
    const bool  commandLineTargetsPmeOnGpu = (pmeOptionArgument == "gpu");
    if (commandLineTargetsPmeOnGpu)
    {
        messages.appendIf(getCompatibleDevices(s_hwinfo->deviceInfoList).empty(),
                          "it targets GPU execution, but no compatible devices were detected");

        if (!commandLineTargetsPmeOnlyRanks && numRanks > 1)
        {
            const bool pmeDecompositionSupported = GMX_USE_cuFFTMp || GMX_USE_Heffte;
            messages.appendIf(!pmeDecompositionSupported,
                              "it targets PME decomposition, but that is not supported");
            if (pmeDecompositionSupported)
            {
                const bool pmeDecompositionActive = (getenv("GMX_GPU_PME_DECOMPOSITION") != nullptr);
                messages.appendIf(!pmeDecompositionActive,
                                  "it targets PME decomposition, but that is not enabled");
                GpuAwareMpiStatus gpuAwareMpiStatus = s_hwinfo->minGpuAwareMpiStatus;
                const bool        gpuAwareMpiActive = gpuAwareMpiStatus == GpuAwareMpiStatus::Forced
                                               || gpuAwareMpiStatus == GpuAwareMpiStatus::Supported;
                messages.appendIf(!gpuAwareMpiActive,
                                  "it targets PME decomposition, which requires GPU-aware MPI, but "
                                  "that is not detected");
            }
        }

        std::optional<std::string_view> pmeFftOptionArgument = commandLine.argumentOf("-pmefft");
        const bool                      commandLineTargetsPmeFftOnGpu =
                !pmeFftOptionArgument.has_value() || pmeFftOptionArgument.value() == "gpu";

        static constexpr bool sc_gpuBuildSyclWithoutGpuFft =
                // NOLINTNEXTLINE(misc-redundant-expression)
                (GMX_GPU_SYCL != 0) && (GMX_GPU_FFT_MKL == 0) && (GMX_GPU_FFT_ROCFFT == 0)
                && (GMX_GPU_FFT_VKFFT == 0) && (GMX_GPU_FFT_BBFFT == 0)
                && (GMX_GPU_FFT_ONEMKL == 0); // NOLINT(misc-redundant-expression)
        messages.appendIf(commandLineTargetsPmeFftOnGpu && sc_gpuBuildSyclWithoutGpuFft,
                          "it targets GPU execution of FFT work, which is not supported in the "
                          "current build");

        std::string errorMessage;
        messages.appendIf(!pme_gpu_supports_build(&errorMessage), errorMessage);
        // A check on whether the .tpr is supported for PME on GPUs is
        // not needed, because it is supported by design.
    }
    return messages;
}

TEST_P(PmeTest, Runs)
{
    auto [pmeTestFlavor, mdrunCommandLine] = GetParam();
    CommandLine commandLine(splitString(mdrunCommandLine));

    // Run mdrun on the tpr file that was built in SetUpTestSuite()
    runner_.tprFileName_ = s_tprFileNames[pmeTestFlavor];

    const bool thisRankChecks = (gmx_node_rank() == 0);
    // Some indentation preserved only for reviewer convenience
    {
        // When using PME tuning on a short mdrun, nstlist needs to be set very short
        // so that the tuning might do something while the test is running.
        const bool usePmeTuning = commandLine.contains("-tunepme");
        if (usePmeTuning)
        {
            commandLine.addOption("-nstlist", 1); // a new grid every step
        }

        ASSERT_TRUE(commandLine.argumentOf("-npme").has_value())
                << "-npme argument is required for this test";
        ASSERT_TRUE(commandLine.argumentOf("-pme").has_value())
                << "-pme argument is required for this test";
        MessageStringCollector skipMessages = getSkipMessagesIfNecessary(commandLine);
        if (!skipMessages.isEmpty())
        {
            GTEST_SKIP() << skipMessages.toString();
        }

        ASSERT_EQ(0, runner_.callMdrun(commandLine));

        if (thisRankChecks)
        {
            // Check the contents of the edr file. Only the main
            // rank should do this I/O intensive operation
            checkEnergies(usePmeTuning);
        }
    }
}

void PmeTest::checkEnergies(const bool usePmeTuning) const
{
    // Some indentation preserved only for reviewer convenience
    {
        {
            TestReferenceData    refData(makeRefDataFileName());
            TestReferenceChecker rootChecker(refData.rootChecker());

            auto energyReader = openEnergyFileToReadTerms(
                    runner_.edrFileName_, { "Coul. recip.", "Total Energy", "Kinetic En." });
            auto conservedChecker  = rootChecker.checkCompound("Energy", "Conserved");
            auto reciprocalChecker = rootChecker.checkCompound("Energy", "Reciprocal");
            // PME tuning causes differing grids and differing
            // reciprocal energy, so we don't check against the same
            // reciprocal energy computed on the CPU
            if (usePmeTuning)
            {
                reciprocalChecker.disableUnusedEntriesCheck();
            }
            bool firstIteration = true;
            while (energyReader->readNextFrame())
            {
                const EnergyFrame& frame            = energyReader->frame();
                const std::string  stepName         = frame.frameName();
                const real         conservedEnergy  = frame.at("Total Energy");
                const real         reciprocalEnergy = frame.at("Coul. recip.");
                if (firstIteration)
                {
                    // use first step values as references for tolerance
                    const real startingKineticEnergy = frame.at("Kinetic En.");
                    const auto conservedTolerance =
                            relativeToleranceAsFloatingPoint(startingKineticEnergy, 2e-5);
                    const auto reciprocalTolerance =
                            relativeToleranceAsFloatingPoint(reciprocalEnergy, 3e-5);
                    reciprocalChecker.setDefaultTolerance(reciprocalTolerance);
                    conservedChecker.setDefaultTolerance(conservedTolerance);
                    firstIteration = false;
                }
                conservedChecker.checkReal(conservedEnergy, stepName.c_str());
                // When not using PME tuning, the reciprocal energy is
                // reproducible enough to check.
                if (!usePmeTuning)
                {
                    reciprocalChecker.checkReal(reciprocalEnergy, stepName.c_str());
                }
            }
        }
    }
}

// To keep test execution time down, we check auto and pme tuning only
// in the Basic case.
//
// Note that some of these cases can only run when there is one MPI
// rank and some require more than one MPI rank. CTest has been
// instructed to run the test binary twice, with respectively one and
// two ranks, so that all tests that can run do run. The test binaries
// consider the hardware, build configuration, and rank count and skip
// those tests that they cannot run.
const auto c_reproducesEnergies = ::testing::ValuesIn(std::vector<PmeTestParameters>{ {
        // Here are all tests without a PME-only rank. These can
        // always run with a single rank, but can only run with two
        // ranks when not targeting GPUs.
        // Note that an -npme argument is required.
        { PmeTestFlavor::Basic, "-notunepme -npme 0 -pme cpu" },
        { PmeTestFlavor::Basic, "-notunepme -npme 0 -pme auto" },
        { PmeTestFlavor::Basic, "-notunepme -npme 0 -pme gpu -pmefft cpu" },
        { PmeTestFlavor::Basic, "-notunepme -npme 0 -pme gpu -pmefft gpu" },
        { PmeTestFlavor::Basic, "-notunepme -npme 0 -pme gpu -pmefft auto" },
        { PmeTestFlavor::WithWalls, "-notunepme -npme 0 -pme cpu" },
        { PmeTestFlavor::WithWalls, "-notunepme -npme 0 -pme gpu -pmefft cpu" },
        { PmeTestFlavor::WithWalls, "-notunepme -npme 0 -pme gpu -pmefft gpu" },
        // Here are all tests with a PME-only rank, which requires
        // more than one total rank
        { PmeTestFlavor::Basic, "-notunepme -npme 1 -pme cpu" },
        { PmeTestFlavor::Basic, "-notunepme -npme 1 -pme auto" },
        { PmeTestFlavor::Basic, "-notunepme -npme 1 -pme gpu -pmefft cpu" },
        { PmeTestFlavor::Basic, "-notunepme -npme 1 -pme gpu -pmefft gpu" },
        { PmeTestFlavor::Basic, "-notunepme -npme 1 -pme gpu -pmefft auto" },
        { PmeTestFlavor::WithWalls, "-notunepme -npme 1 -pme cpu" },
        { PmeTestFlavor::WithWalls, "-notunepme -npme 1 -pme gpu -pmefft cpu" },
        { PmeTestFlavor::WithWalls, "-notunepme -npme 1 -pme gpu -pmefft gpu" },
        // All tests with PME tuning here
        { PmeTestFlavor::Basic, "-tunepme -npme 0 -pme cpu" },
        { PmeTestFlavor::Basic, "-tunepme -npme 0 -pme gpu -pmefft cpu" },
        { PmeTestFlavor::Basic, "-tunepme -npme 0 -pme gpu -pmefft gpu" },
} });

INSTANTIATE_TEST_SUITE_P(ReproducesEnergies, PmeTest, c_reproducesEnergies, nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
