/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Test for MD with orientation restraints
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

const std::string g_mdpContents = R"(
dt            = 0.002
nsteps        = 10
tcoupl        = V-rescale
tc-grps       = System
tau-t         = 0.5
ref-t         = 300
constraints   = h-bonds
cutoff-scheme = Verlet
orire         = Yes
orire-fitgrp  = backbone)";

class OriresTest : public MdrunTestFixture
{
};

/* Check whether the orires function works. */
TEST_F(OriresTest, OriresCanRun)
{
    GMX_MPI_TEST(RequireRankCount<1>);

    runner_.useTopGroAndNdxFromDatabase("orires_1lvz");
    runner_.useStringAsMdpFile(g_mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    // Do an mdrun with orientation restraints enabled
    ASSERT_EQ(0, runner_.callMdrun());
}

//! A simple test fixture for a multi-simulation with one rank per simulation
class SimpleMultiSimTest : public MdrunTestFixture
{
public:
    SimpleMultiSimTest()
    {
        // Only single-rank simulations are supported
        const int numRanksPerSimulation = 1;

        {
            // How many MPI ranks exist?
            const int commSize = getNumberOfTestMpiRanks();
            const int commRank = gmx_node_rank();
            // Set up the multi-sim
            ensembleSize_ = commSize / numRanksPerSimulation;
            simulationId_ = commRank / numRanksPerSimulation;
            isMainRank_   = commRank % numRanksPerSimulation == 0;
        }

        // Multi-sim requires organizing creating subdirectories
        // and tpr files inside them.
        const std::filesystem::path& originalTempDirectory = fileManager_.getOutputTempDirectory();

        // Prepare the mdrun caller and make the subdirectories
        mdrunCaller_.append("mdrun");
        mdrunCaller_.addOption("-multidir");
        for (int i = 0; i < ensembleSize_; ++i)
        {
            std::filesystem::path newTempDirectory = originalTempDirectory;
            newTempDirectory.append(formatString("sim_%d", i));
            mdrunCaller_.append(newTempDirectory.string());
            if (i == simulationId_)
            {
                if (isMainRank_)
                {
                    // Only one rank per simulation should actually create the directory!
                    std::filesystem::create_directory(newTempDirectory);
                }
#if GMX_LIB_MPI
                // Make sure the directory has been made before other ranks try to use it as the temp directory.
                MPI_Barrier(MdrunTestFixtureBase::s_communicator);
#endif
                // Use the new directory for output.
                fileManager_.setOutputTempDirectory(newTempDirectory);
                // Replace the SimulationRunner with one using the new
                // output directory. This makes sure grompp and mdrun
                // write their output where it is expected for a
                // multi-sim.
                runner_ = SimulationRunner(&fileManager_);
            }
        }
    }

    //! Number of simulations in the ensemble
    int ensembleSize_;
    //! ID of this simulation within the set
    int simulationId_;
    //!  Whether this rank will be the main rank of a simulation
    bool isMainRank_ = false;
    //! The helper object that will call mdrun -multidir
    CommandLine mdrunCaller_;
};

using OriresEnsembleTest = SimpleMultiSimTest;

TEST_F(OriresEnsembleTest, OriresCanAverageEnsembles)
{
    GMX_MPI_TEST(RequireRankCount<2>);

    if (!GMX_LIB_MPI)
    {
        GTEST_SKIP() << "Library MPI build configuration required for ensemble restraints";
    }
    if (ensembleSize_ <= 1)
    {
        GTEST_SKIP() << "Must have at least two simulations for ensemble restraints";
    }

    // Only one rank per simulation creates mdp and tpr files
    if (isMainRank_)
    {
        runner_.useTopGroAndNdxFromDatabase("orires_1lvz");
        runner_.useStringAsMdpFile(g_mdpContents);
        EXPECT_EQ(0, runner_.callGromppOnThisRank());
    }

    // Do a multi-sim mdrun with orientation restraints enabled
    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller_));
}

} // namespace test
} // namespace gmx
