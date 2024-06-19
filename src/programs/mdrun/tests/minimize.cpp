/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests for the energy minimization functionality.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <cstdio>

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Test fixture base for energy minimizaiton
 *
 * This test ensures mdrun can run an energy minimization, reaching
 * a reproducible final energy.
 *
 * The choices for tolerance are arbitrary but sufficient. */
class EnergyMinimizationTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string>>
{
};

/*! \brief Database of empirical tolerances for EM integrators on the various systems. */
std::unordered_map<std::string, FloatingPointTolerance> potentialEnergyToleranceForSystem_g = {
    { { "argon12", relativeToleranceAsPrecisionDependentUlp(-1, 10, 200) },
      { "tip3p5", relativeToleranceAsPrecisionDependentUlp(-50, 200, 3800) },
      { "glycine_vacuo", relativeToleranceAsPrecisionDependentUlp(1000, 100, 100) },
      { "alanine_vsite_vacuo", relativeToleranceAsPrecisionDependentUlp(-160, 150, 400) },
      { "glycine_no_constraints_vacuo", relativeToleranceAsPrecisionDependentUlp(2000, 100, 100) } }
};

TEST_P(EnergyMinimizationTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto minimizer      = std::get<1>(params);
    SCOPED_TRACE(formatString("Comparing '%s' energy minimization for simulation '%s'",
                              minimizer.c_str(),
                              simulationName.c_str()));

    // TODO At some point we should also test PME-only ranks.
    int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        fprintf(stdout,
                "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName.c_str(),
                numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
        return;
    }

    auto mdpFieldValues =
            prepareMdpFieldValues(simulationName.c_str(), minimizer.c_str(), "no", "no");
    mdpFieldValues["nsteps"] = "4";

    int maxWarningsTolerated = (minimizer == "l-bfgs") ? 1 : 0;
    // prepare the .tpr file
    {
        // TODO evolve grompp to report the number of warnings issued, so
        // tests always expect the right number.
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-maxwarn", maxWarningsTolerated);
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }

    // do mdrun, preparing to check the energies later
    runner_.edrFileName_ = fileManager_.getTemporaryFilePath("minimize.edr").string();
    {
        CommandLine mdrunCaller;
        mdrunCaller.append("mdrun");
        if (minimizer == "l-bfgs" && getNumberOfTestMpiRanks() > 1)
        {
            // Ideally we would use this death test, but it is not
            // stable enough in CI, so we just skip it.
            // EXPECT_DEATH_IF_SUPPORTED(runner_.callMdrun(mdrunCaller),
            //                          "L-BFGS minimization only supports a single rank");
            return;
        }
        else
        {
            ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
        }
    }

    EnergyTermsToCompare energyTermsToCompare{ {
            { interaction_function[F_EPOT].longname, potentialEnergyToleranceForSystem_g.at(simulationName) },
    } };

    TestReferenceData refData;
    auto              checker = refData.rootChecker()
                           .checkCompound("Simulation", simulationName)
                           .checkCompound("Minimizer", minimizer);
    checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
}

//! Containers of systems and integrators to test.
//! \{
std::vector<std::string> unconstrainedSystemsToTest_g = { "argon12",
                                                          "glycine_no_constraints_vacuo" };
std::vector<std::string> minimizersToTest_g           = { "steep", "cg", "l-bfgs" };

std::vector<std::string> constrainedSystemsToTest_g        = { "tip3p5",
                                                        "glycine_vacuo",
                                                        "alanine_vsite_vacuo" };
std::vector<std::string> minimizersToTestWithConstraints_g = { "steep", "cg" };
//! \}

// The time for OpenCL kernel compilation means these tests might time
// out. If that proves to be a problem, these can be disabled for
// OpenCL builds. However, once that compilation is cached for the
// lifetime of the whole test binary process, these tests should run in
// such configurations.
INSTANTIATE_TEST_SUITE_P(MinimizersWorkWithConstraints,
                         EnergyMinimizationTest,
                         ::testing::Combine(::testing::ValuesIn(constrainedSystemsToTest_g),
                                            ::testing::ValuesIn(minimizersToTestWithConstraints_g)));
INSTANTIATE_TEST_SUITE_P(MinimizersWork,
                         EnergyMinimizationTest,
                         ::testing::Combine(::testing::ValuesIn(unconstrainedSystemsToTest_g),
                                            ::testing::ValuesIn(minimizersToTest_g)));

} // namespace
} // namespace test
} // namespace gmx
