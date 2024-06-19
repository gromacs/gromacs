/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests for the mdrun multi-simulation functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "multisimtest.h"

#include "config.h"

#include <cmath>

#include <algorithm>
#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"

#include "moduletest.h"
#include "terminationhelper.h"

namespace gmx
{
namespace test
{

MultiSimTest::MultiSimTest() :
    size_(gmx_node_num()),
    rank_(gmx_node_rank()),
    numRanksPerSimulation_(std::get<0>(GetParam())),
    simulationNumber_(rank_ / numRanksPerSimulation_),
    mdrunCaller_(new CommandLine)

{
    // Zero or less ranks doesn't make sense
    GMX_RELEASE_ASSERT(numRanksPerSimulation_ > 0, "Invalid number of ranks per simulation.");

    const char* directoryNameFormat = "sim_%d";

    // Modify the file manager to have a temporary directory unique to
    // each simulation. No need to have a mutex on this, nobody else
    // can access the fileManager_ yet because we only just
    // constructed it.
    const std::filesystem::path& originalTempDirectory = fileManager_.getOutputTempDirectory();
    const std::filesystem::path  newTempDirectory =
            std::filesystem::path(originalTempDirectory).append(formatString(directoryNameFormat, simulationNumber_));
    if (rank_ % numRanksPerSimulation_ == 0)
    {
        // Only one rank per simulation creates directory
        std::filesystem::create_directory(newTempDirectory);
    }
#if GMX_LIB_MPI
    // Make sure directories got created.
    MPI_Barrier(MdrunTestFixtureBase::s_communicator);
#endif
    fileManager_.setOutputTempDirectory(newTempDirectory);

    mdrunCaller_->append("mdrun");
    mdrunCaller_->addOption("-multidir");
    for (int i = 0; i < size_ / numRanksPerSimulation_; ++i)
    {
        mdrunCaller_->append(std::filesystem::path(originalTempDirectory)
                                     .append(formatString(directoryNameFormat, i))
                                     .string());
    }
}

bool MultiSimTest::mpiSetupValid() const
{
    // Single simulation case is not implemented in multi-sim
    const bool haveAtLeastTwoSimulations = ((size_ / numRanksPerSimulation_) >= 2);
    // Mdrun will throw error if simulations don't have identical number of ranks
    const bool simulationsHaveIdenticalRankNumber = ((size_ % numRanksPerSimulation_) == 0);

    return (haveAtLeastTwoSimulations && simulationsHaveIdenticalRankNumber);
}

void MultiSimTest::organizeMdpFile(SimulationRunner*    runner,
                                   IntegrationAlgorithm integrator,
                                   TemperatureCoupling  tcoupl,
                                   PressureCoupling     pcoupl,
                                   int                  numSteps,
                                   bool                 doRegression) const
{
    GMX_RELEASE_ASSERT(mpiSetupValid(), "Creating the mdp file without valid MPI setup is useless.");
    const real  baseTemperature = 298;
    const real  basePressure    = 1;
    std::string mdpFileContents = formatString(
            "integrator = %s\n"
            "tcoupl = %s\n"
            "nsttcouple = 10\n"
            "pcoupl = %s\n"
            "nstpcouple = 10\n"
            "nsteps = %d\n"
            "nstlog = 1\n"
            "nstcalcenergy = 1\n"
            "tc-grps = System\n"
            "tau-t = 1\n"
            "ref-t = %f\n"
            // pressure coupling (if active)
            "tau-p = 2\n"
            "ref-p = %f\n"
            "compressibility = 4.5e-5\n"
            // velocity generation
            "gen-vel = yes\n"
            "gen-temp = %f\n"
            "gen-seed = %d\n"
            // v-rescale and c-rescale use ld-seed
            "ld-seed = %d\n"
            // Two systems are used: spc2    non-interacting also at cutoff 1nm
            //                       tip3p5  has box length of 1.86
            "rcoulomb = 0.7\n"
            "rvdw = 0.7\n"
            // Trajectory output if required
            "nstxout = %d\n"
            "nstvout = %d\n"
            "nstfout = %d\n"
            "nstenergy = %d\n",
            enumValueToString(integrator),
            enumValueToString(tcoupl),
            enumValueToString(pcoupl),
            numSteps,
            baseTemperature + 0.0001 * rank_,
            basePressure * std::pow(1.01, rank_),
            /* Set things up so that the initial KE decreases with
               increasing replica number, so that the (identical)
               starting PE decreases on the first step more for the
               replicas with higher number, which will tend to force
               replica exchange to occur. */
            std::max(baseTemperature - 10 * rank_, real(0)),
            // If we do regression, we need reproducible velocity
            // generation, which can be different per simulation
            (doRegression ? 671324 + simulationNumber_ : -1),
            // If we do regression, we need reproducible temperature and
            // pressure coupling, which can be different per simulation
            (doRegression ? 51203 + simulationNumber_ : -1),
            // If we do regression, write one intermediate point
            (doRegression ? int(numSteps / 2) : 0),
            (doRegression ? int(numSteps / 2) : 0),
            (doRegression ? int(numSteps / 2) : 0),
            // If we do regression, print energies every step so
            // we're sure to catch the replica exchange steps
            (doRegression ? 1 : 1000));
    runner->useStringAsMdpFile(mdpFileContents);
}

void MultiSimTest::runGrompp(SimulationRunner* runner, int numSteps, bool doRegression, int maxWarnings) const
{
    // Call grompp once per simulation
    if (rank_ % numRanksPerSimulation_ == 0)
    {
        const auto& simulator = std::get<1>(GetParam());
        const auto& tcoupl    = std::get<2>(GetParam());
        const auto& pcoupl    = std::get<3>(GetParam());
        int         maxWarn   = maxWarnings;
        if (pcoupl == PressureCoupling::Berendsen)
        {
            maxWarn++;
        }
        if (tcoupl == TemperatureCoupling::Berendsen)
        {
            maxWarn++;
        }
        organizeMdpFile(runner, simulator, tcoupl, pcoupl, numSteps, doRegression);
        runner->setMaxWarn(maxWarn);
        CommandLine caller;
        EXPECT_EQ(0, runner->callGromppOnThisRank(caller));
    }

#if GMX_LIB_MPI
    // Make sure simulation mains have written the .tpr file before other ranks try to read it.
    MPI_Barrier(MdrunTestFixtureBase::s_communicator);
#endif
}

void MultiSimTest::runExitsNormallyTest()
{
    if (!mpiSetupValid())
    {
        // Can't test multi-sim without multiple simulations
        return;
    }

    SimulationRunner runner(&fileManager_);
    runner.useTopGroAndNdxFromDatabase("spc2");

    runGrompp(&runner);

    ASSERT_EQ(0, runner.callMdrun(*mdrunCaller_));
}

void MultiSimTest::runMaxhTest()
{
    if (!mpiSetupValid())
    {
        // Can't test multi-sim without multiple simulations
        return;
    }

    SimulationRunner runner(&fileManager_);
    runner.useTopGroAndNdxFromDatabase("spc2");

    TerminationHelper helper(mdrunCaller_.get(), &runner);
    // Make sure -maxh has a chance to propagate
    int numSteps = 100;
    runGrompp(&runner, numSteps);

    helper.runFirstMdrun(runner.cptOutputFileName_);
    helper.runSecondMdrun();
}

} // namespace test
} // namespace gmx
