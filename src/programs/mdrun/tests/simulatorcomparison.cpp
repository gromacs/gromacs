/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Helper classes for tests that compare the results of equivalent
 * simulation runs. Currently used for the rerun and the simulator
 * tests
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "testutils/mpitest.h"
#include "testutils/setenv.h"
#include "testutils/simulationdatabase.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "mdruncomparison.h"
#include "moduletest.h"
#include "trajectorycomparison.h"
#include "trajectoryreader.h"

namespace gmx
{
namespace test
{
namespace
{

//! Run grompp and mdrun for both sets of mdp field values
template<bool doEnvironmentVariable, bool doRerun>
void executeSimulatorComparisonTestImpl(TestFileManager*            fileManager,
                                        SimulationRunner*           runner,
                                        const std::string&          simulationName,
                                        int                         maxWarningsTolerated,
                                        const MdpFieldValues&       mdpFieldValues,
                                        const EnergyTermsToCompare& energyTermsToCompare,
                                        const TrajectoryComparison& trajectoryComparison,
                                        const std::string&          environmentVariable)
{
    // TODO At some point we should also test PME-only ranks.
    int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        fprintf(stdout,
                "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName.c_str(), numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
        return;
    }

    auto simulator1TrajectoryFileName = fileManager->getTemporaryFilePath("sim1.trr");
    auto simulator1EdrFileName        = fileManager->getTemporaryFilePath("sim1.edr");
    auto simulator2TrajectoryFileName = fileManager->getTemporaryFilePath("sim2.trr");
    auto simulator2EdrFileName        = fileManager->getTemporaryFilePath("sim2.edr");
    auto simulatorTprFileName         = fileManager->getTemporaryFilePath("sim.tpr");

    // prepare the .tpr file
    {
        // TODO evolve grompp to report the number of warnings issued, so
        // tests always expect the right number.
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-maxwarn", maxWarningsTolerated);
        runner->tprFileName_ = simulatorTprFileName;
        runner->useTopGroAndNdxFromDatabase(simulationName);
        runner->useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner->callGrompp(caller));
    }

    char* environmentVariableBackup = nullptr;
    if (doEnvironmentVariable)
    {
        // save state of environment variable
        environmentVariableBackup = getenv(environmentVariable.c_str());
    }

    // do the first mdrun
    {
        runner->fullPrecisionTrajectoryFileName_ = simulator1TrajectoryFileName;
        runner->edrFileName_                     = simulator1EdrFileName;
        runner->tprFileName_                     = simulatorTprFileName;
        CommandLine simulator1Caller;
        simulator1Caller.append("mdrun");
        if (doEnvironmentVariable)
        {
            // unset environment variable
            gmxUnsetenv(environmentVariable.c_str());
        }
        ASSERT_EQ(0, runner->callMdrun(simulator1Caller));
    }

    // do the second mdrun
    {
        runner->fullPrecisionTrajectoryFileName_ = simulator2TrajectoryFileName;
        runner->edrFileName_                     = simulator2EdrFileName;
        runner->tprFileName_                     = simulatorTprFileName;
        CommandLine simulator2Caller;
        simulator2Caller.append("mdrun");
        if (doEnvironmentVariable)
        {
            // set environment variable
            gmxSetenv(environmentVariable.c_str(), "ON", true);
        }
        if (doRerun)
        {
            simulator2Caller.addOption("-rerun", simulator1TrajectoryFileName);
        }
        ASSERT_EQ(0, runner->callMdrun(simulator2Caller));
    }

    if (doEnvironmentVariable)
    {
        if (environmentVariableBackup != nullptr)
        {
            // set environment variable
            gmxSetenv(environmentVariable.c_str(), environmentVariableBackup, true);
        }
        else
        {
            // unset environment variable
            gmxUnsetenv(environmentVariable.c_str());
        }
    }

    // Build the functor that will compare energy frames on the chosen
    // energy terms.
    EnergyComparison energyComparison(energyTermsToCompare);

    // Build the manager that will present matching pairs of frames to compare.
    //
    // TODO Here is an unnecessary copy of keys (ie. the energy term
    // names), for convenience. In the future, use a range.
    auto                                namesOfEnergiesToMatch = energyComparison.getEnergyNames();
    FramePairManager<EnergyFrameReader> energyManager(
            openEnergyFileToReadTerms(simulator1EdrFileName, namesOfEnergiesToMatch),
            openEnergyFileToReadTerms(simulator2EdrFileName, namesOfEnergiesToMatch));
    // Compare the energy frames.
    energyManager.compareAllFramePairs<EnergyFrame>(energyComparison);

    // Build the manager that will present matching pairs of frames to compare
    FramePairManager<TrajectoryFrameReader> trajectoryManager(
            std::make_unique<TrajectoryFrameReader>(simulator1TrajectoryFileName),
            std::make_unique<TrajectoryFrameReader>(simulator2TrajectoryFileName));
    // Compare the trajectory frames.
    trajectoryManager.compareAllFramePairs<TrajectoryFrame>(trajectoryComparison);
}
} // namespace

template<typename... Args>
void executeSimulatorComparisonTest(const std::string& environmentVariable, Args&&... args)
{
    executeSimulatorComparisonTestImpl<true, false>(std::forward<Args>(args)..., environmentVariable);
}

template<typename... Args>
void executeRerunTest(Args&&... args)
{
    executeSimulatorComparisonTestImpl<false, true>(std::forward<Args>(args)..., "");
}

} // namespace test
} // namespace gmx
