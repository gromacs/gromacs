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

#include "config.h"

#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/mdtypes/md_enums.h"

#include "programs/mdrun/tests/moduletest.h"

#include "multisimtest.h"

namespace gmx
{
namespace test
{

/* This test ensures mdrun can run multi-simulations.  It runs one
 * simulation per MPI rank.
 *
 * TODO Preferably, we could test that mdrun correctly refuses to run
 * multi-simulation unless compiled with real MPI with more than one
 * rank available. However, if we just call mdrun blindly, those cases
 * trigger an error that is currently fatal to mdrun and also to the
 * test binary. So, in the meantime we must not test those cases. If
 * there is no MPI, we disable the test, so that there is a reminder
 * that it is disabled. There's no elegant way to conditionally
 * disable a test at run time, so currently there is no feedback if
 * only one rank is available. However, the test harness knows to run
 * this test with more than one rank. */
TEST_P(MultiSimTest, ExitsNormally)
{
    runExitsNormallyTest();
}

TEST_P(MultiSimTest, ExitsNormallyWithDifferentNumbersOfStepsPerSimulation)
{
    if (!mpiSetupValid())
    {
        // MPI setup is not suitable for multi-sim
        return;
    }
    SimulationRunner runner(&fileManager_);
    runner.useTopGroAndNdxFromDatabase("spc2");

    // Do some different small numbers of steps in each simulation
    int numSteps = simulationNumber_ % 4;
    runGrompp(&runner, numSteps);

    ASSERT_EQ(0, runner.callMdrun(*mdrunCaller_));
}

/* Note, not all preprocessor implementations nest macro expansions
   the same way / at all, if we would try to duplicate less code. */
#if GMX_LIB_MPI
INSTANTIATE_TEST_SUITE_P(InNvt,
                         MultiSimTest,
                         ::testing::Combine(::testing::Values(NumRanksPerSimulation(1),
                                                              NumRanksPerSimulation(2)),
                                            ::testing::Values(IntegrationAlgorithm::MD),
                                            ::testing::Values(TemperatureCoupling::VRescale),
                                            ::testing::Values(PressureCoupling::No)));
#else
// Test needs real MPI to run
INSTANTIATE_TEST_SUITE_P(DISABLED_InNvt,
                         MultiSimTest,
                         ::testing::Combine(::testing::Values(NumRanksPerSimulation(1),
                                                              NumRanksPerSimulation(2)),
                                            ::testing::Values(IntegrationAlgorithm::MD),
                                            ::testing::Values(TemperatureCoupling::VRescale),
                                            ::testing::Values(PressureCoupling::No)));
#endif

//! Convenience typedef
typedef MultiSimTest MultiSimTerminationTest;

TEST_P(MultiSimTerminationTest, WritesCheckpointAfterMaxhTerminationAndThenRestarts)
{
    runMaxhTest();
}

INSTANTIATE_TEST_SUITE_P(InNvt,
                         MultiSimTerminationTest,
                         ::testing::Combine(::testing::Values(NumRanksPerSimulation(1),
                                                              NumRanksPerSimulation(2)),
                                            ::testing::Values(IntegrationAlgorithm::MD),
                                            ::testing::Values(TemperatureCoupling::VRescale),
                                            ::testing::Values(PressureCoupling::No)));

} // namespace test
} // namespace gmx
