/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2020, by the GROMACS development team.
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * This implements basic PME sanity tests.
 * It runs the input system with PME for several steps (on CPU and GPU, if available),
 * and checks the reciprocal and conserved energies.
 * As part of mdrun-test, this will always run single rank PME simulation.
 * As part of mdrun-mpi-test, this will run same as above when a single rank is requested,
 * or a simulation with a single separate PME rank ("-npme 1") when multiple ranks are requested.
 * \todo Extend and generalize this for more multi-rank tests (-npme 0, -npme 2, etc).
 * \todo Implement death tests (e.g. for PME GPU decomposition).
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <map>
#include <mutex>
#include <string>
#include <vector>

#include <gtest/gtest-spi.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/refdata.h"

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

/*! \brief A basic PME runner
 *
 * \todo Consider also using GpuTest class. */
class PmeTest : public MdrunTestFixture
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
    //! Convenience typedef
    using RunModesList = std::map<std::string, std::vector<const char*>>;
    //! Runs the test with the given inputs
    void runTest(const RunModesList& runModes, PmeTestFlavor pmeTestFlavor);
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
        runner.tprFileName_ = s_testFileManager->getTemporaryFilePath(tprFileNameSuffix);
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
        messages.appendIf(!commandLineTargetsPmeOnlyRanks && numRanks > 1,
                          "it targets PME decomposition, but that is not supported");

        std::optional<std::string_view> pmeFftOptionArgument = commandLine.argumentOf("-pmefft");
        const bool                      commandLineTargetsPmeFftOnGpu =
                !pmeFftOptionArgument.has_value() || pmeFftOptionArgument.value() == "gpu";
        messages.appendIf(commandLineTargetsPmeFftOnGpu && GMX_GPU_SYCL, // Issues #4219, #4274
                          "it targets GPU execution of FFT work, which is not supported with DPC++ "
                          "or hipSYCL");

        std::string errorMessage;
        messages.appendIf(!pme_gpu_supports_build(&errorMessage), errorMessage);
        messages.appendIf(!pme_gpu_supports_hardware(*s_hwinfo, &errorMessage), errorMessage);
        // A check on whether the .tpr is supported for PME on GPUs is
        // not needed, because it is supported by design.
    }
    return messages;
}

void PmeTest::runTest(const RunModesList& runModes, const PmeTestFlavor pmeTestFlavor)
{
    const std::string inputFile = "spc-and-methanol";
    runner_.useTopGroAndNdxFromDatabase(inputFile);

    // With single rank we can and will always test PP+PME as part of mdrun-test.
    // With multiple ranks we can additionally test a single PME-only rank within mdrun-mpi-test.
    const bool parallelRun    = (getNumberOfTestMpiRanks() > 1);
    const bool useSeparatePme = parallelRun;

    // Run mdrun on the tpr file that was built in SetUpTestSuite()
    runner_.tprFileName_ = s_tprFileNames[pmeTestFlavor];

    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());
    const bool           thisRankChecks = (gmx_node_rank() == 0);
    if (!thisRankChecks)
    {
        EXPECT_NONFATAL_FAILURE(rootChecker.checkUnusedEntries(), ""); // skip checks on other ranks
    }

    for (const auto& mode : runModes)
    {
        SCOPED_TRACE("mdrun " + joinStrings(mode.second, " "));
        runner_.edrFileName_ =
                fileManager_.getTemporaryFilePath(inputFile + "_" + mode.first + ".edr");

        CommandLine commandLine(mode.second);

        const bool usePmeTuning = (mode.first.find("Tune") != std::string::npos);
        if (usePmeTuning)
        {
            commandLine.append("-tunepme");
            commandLine.addOption("-nstlist", 1); // a new grid every step
        }
        else
        {
            commandLine.append("-notunepme"); // for reciprocal energy reproducibility
        }
        if (useSeparatePme)
        {
            commandLine.addOption("-npme", 1);
        }
        else
        {
            commandLine.addOption("-npme", 0);
        }

        ASSERT_TRUE(commandLine.argumentOf("-npme").has_value())
                << "-npme argument is required for this test";
        ASSERT_TRUE(commandLine.argumentOf("-pme").has_value())
                << "-pme argument is required for this test";
        MessageStringCollector skipMessages = getSkipMessagesIfNecessary(commandLine);
        if (!skipMessages.isEmpty())
        {
            // For now, just write these to stdout. Later, we can use
            // GTEST_SKIP().
            fputs(skipMessages.toString().c_str(), stdout);
            continue;
        }

        ASSERT_EQ(0, runner_.callMdrun(commandLine));

        if (thisRankChecks)
        {
            auto energyReader = openEnergyFileToReadTerms(
                    runner_.edrFileName_, { "Coul. recip.", "Total Energy", "Kinetic En." });
            auto conservedChecker  = rootChecker.checkCompound("Energy", "Conserved");
            auto reciprocalChecker = rootChecker.checkCompound("Energy", "Reciprocal");
            bool firstIteration    = true;
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
                if (!usePmeTuning) // with PME tuning come differing grids and differing reciprocal energy
                {
                    reciprocalChecker.checkReal(reciprocalEnergy, stepName.c_str());
                }
            }
        }
    }
}

TEST_F(PmeTest, ReproducesEnergies)
{
    // TODO test all proper/improper combinations in more thorough way?
    RunModesList runModes;
    runModes["PmeOnCpu"]         = { "-pme", "cpu" };
    runModes["PmeAuto"]          = { "-pme", "auto" };
    runModes["PmeOnGpuFftOnCpu"] = { "-pme", "gpu", "-pmefft", "cpu" };
    runModes["PmeOnGpuFftOnGpu"] = { "-pme", "gpu", "-pmefft", "gpu" };
    runModes["PmeOnGpuFftAuto"]  = { "-pme", "gpu", "-pmefft", "auto" };
    // same manual modes but marked for PME tuning
    runModes["PmeOnCpuTune"]         = { "-pme", "cpu" };
    runModes["PmeOnGpuFftOnCpuTune"] = { "-pme", "gpu", "-pmefft", "cpu" };
    runModes["PmeOnGpuFftOnGpuTune"] = { "-pme", "gpu", "-pmefft", "gpu" };

    runTest(runModes, PmeTestFlavor::Basic);
}

TEST_F(PmeTest, ScalesTheBoxWithWalls)
{
    RunModesList runModes;
    runModes["PmeOnCpu"]         = { "-pme", "cpu" };
    runModes["PmeOnGpuFftOnCpu"] = { "-pme", "gpu", "-pmefft", "cpu" };
    runModes["PmeOnGpuFftOnGpu"] = { "-pme", "gpu", "-pmefft", "gpu" };

    runTest(runModes, PmeTestFlavor::WithWalls);
}

} // namespace
} // namespace test
} // namespace gmx
