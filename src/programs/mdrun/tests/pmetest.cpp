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
#include "gromacs/utility/gmxmpi.h"
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

/*! \brief A basic PME runner
 *
 * \todo Consider also using GpuTest class. */
class PmeTest : public MdrunTestFixture
{
public:
    //! Convenience typedef
    using RunModesList = std::map<std::string, std::vector<const char*>>;
    //! Runs the test with the given inputs
    void runTest(const RunModesList& runModes);
};

void PmeTest::runTest(const RunModesList& runModes)
{
    const std::string inputFile = "spc-and-methanol";
    runner_.useTopGroAndNdxFromDatabase(inputFile);

    // With single rank we can and will always test PP+PME as part of mdrun-test.
    // With multiple ranks we can additionally test a single PME-only rank within mdrun-mpi-test.
    const bool parallelRun    = (getNumberOfTestMpiRanks() > 1);
    const bool useSeparatePme = parallelRun;

    EXPECT_EQ(0, runner_.callGrompp());

    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());
    const bool           thisRankChecks = (gmx_node_rank() == 0);
    if (!thisRankChecks)
    {
        EXPECT_NONFATAL_FAILURE(rootChecker.checkUnusedEntries(), ""); // skip checks on other ranks
    }

    auto hardwareInfo_ =
            gmx_detect_hardware(PhysicalNodeCommunicator(MPI_COMM_WORLD, gmx_physicalnode_id_hash()));

    for (const auto& mode : runModes)
    {
        SCOPED_TRACE("mdrun " + joinStrings(mode.second, " "));
        auto modeTargetsGpus = (mode.first.find("Gpu") != std::string::npos);
        if (modeTargetsGpus && getCompatibleDevices(hardwareInfo_->deviceInfoList).empty())
        {
            // This run mode will cause a fatal error from mdrun when
            // it can't find GPUs, which is not something we're trying
            // to test here.
            continue;
        }
        auto modeTargetsPmeOnGpus = (mode.first.find("PmeOnGpu") != std::string::npos);
        if (modeTargetsPmeOnGpus
            && !(pme_gpu_supports_build(nullptr) && pme_gpu_supports_hardware(*hardwareInfo_, nullptr)))
        {
            // This run mode will cause a fatal error from mdrun when
            // it finds an unsuitable device, which is not something
            // we're trying to test here.
            continue;
        }

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
    const int         nsteps     = 20;
    const std::string theMdpFile = formatString(
            "coulombtype     = PME\n"
            "nstcalcenergy   = 1\n"
            "nstenergy       = 1\n"
            "pme-order       = 4\n"
            "nsteps          = %d\n",
            nsteps);

    runner_.useStringAsMdpFile(theMdpFile);

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

    runTest(runModes);
}

TEST_F(PmeTest, ScalesTheBox)
{
    const int         nsteps     = 0;
    const std::string theMdpFile = formatString(
            "coulombtype     = PME\n"
            "nstcalcenergy   = 1\n"
            "nstenergy       = 1\n"
            "pme-order       = 4\n"
            "pbc             = xyz\n"
            "nsteps          = %d\n",
            nsteps);

    runner_.useStringAsMdpFile(theMdpFile);

    RunModesList runModes;
    runModes["PmeOnCpu"]         = { "-pme", "cpu" };
    runModes["PmeOnGpuFftOnCpu"] = { "-pme", "gpu", "-pmefft", "cpu" };
    runModes["PmeOnGpuFftOnGpu"] = { "-pme", "gpu", "-pmefft", "gpu" };

    runTest(runModes);
}

TEST_F(PmeTest, ScalesTheBoxWithWalls)
{
    const int         nsteps     = 0;
    const std::string theMdpFile = formatString(
            "coulombtype     = PME\n"
            "nstcalcenergy   = 1\n"
            "nstenergy       = 1\n"
            "pme-order       = 4\n"
            "pbc             = xy\n"
            "nwall           = 2\n"
            "ewald-geometry  = 3dc\n"
            "wall_atomtype   = CMet H\n"
            "wall_density    = 9 9.0\n"
            "wall-ewald-zfac = 5\n"
            "nsteps          = %d\n",
            nsteps);

    runner_.useStringAsMdpFile(theMdpFile);

    RunModesList runModes;
    runModes["PmeOnCpu"]         = { "-pme", "cpu" };
    runModes["PmeOnGpuFftOnCpu"] = { "-pme", "gpu", "-pmefft", "cpu" };
    runModes["PmeOnGpuFftOnGpu"] = { "-pme", "gpu", "-pmefft", "gpu" };

    runTest(runModes);
}

} // namespace
} // namespace test
} // namespace gmx
