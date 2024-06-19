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
 * Tests for the normal modes functionality.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <cstdio>

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/filenameoption.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/xvgtest.h"

#include "energycomparison.h"
#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Test fixture base for normal mode analysis
 *
 * This test ensures mdrun can run a normal mode analysis, reaching
 * a reproducible eigenvalues following diagonalization.
 *
 * The choices for tolerance are arbitrary but sufficient. */
class NormalModesTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string>>
{
};

TEST_P(NormalModesTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto integrator     = std::get<1>(params);
    SCOPED_TRACE(formatString("Comparing normal modes for '%s'", simulationName.c_str()));

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
            prepareMdpFieldValues(simulationName.c_str(), integrator.c_str(), "no", "no");
    mdpFieldValues["nsteps"]      = "1";
    mdpFieldValues["rcoulomb"]    = "5.6";
    mdpFieldValues["rlist"]       = "5.6";
    mdpFieldValues["rvdw"]        = "5.6";
    mdpFieldValues["constraints"] = "none";
    mdpFieldValues["coulombtype"] = "Cut-off";
    mdpFieldValues["vdwtype"]     = "Cut-off";

    // prepare the .tpr file
    {
        CommandLine caller;
        runner_.useTopG96AndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Do mdrun, preparing to check the normal modes later
    {
        CommandLine mdrunCaller;
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
    }
    // Now run gmx nmeig and check the output
    {
        ASSERT_EQ(0, runner_.callNmeig());
        TestReferenceData refData;
        auto              checker = refData.rootChecker()
                               .checkCompound("System", simulationName)
                               .checkCompound("Integrator", integrator);
        auto settings      = XvgMatchSettings();
        settings.tolerance = relativeToleranceAsFloatingPoint(1.0, 1e-05);
        TextInputFile input("eigenval.xvg");
        checkXvgFile(&input, &checker, settings);
    }
}

// The time for OpenCL kernel compilation means these tests might time
// out. If that proves to be a problem, these can be disabled for
// OpenCL builds. However, once that compilation is cached for the
// lifetime of the whole test binary process, these tests should run in
// such configurations.
#if GMX_DOUBLE

//! Containers of systems and integrators to test.
//! \{
std::vector<std::string> systemsToTest_g     = { "scaled-water",
                                             "villin",
                                             "spc-dimer",
                                             "one-tip5p",
                                             "sw-dimer" };
std::vector<std::string> integratorsToTest_g = { "nm" };

//! \}

INSTANTIATE_TEST_SUITE_P(NormalModesWorks,
                         NormalModesTest,
                         ::testing::Combine(::testing::ValuesIn(systemsToTest_g),
                                            ::testing::ValuesIn(integratorsToTest_g)));
#else
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(NormalModesTest);
#endif
} // namespace
} // namespace test
} // namespace gmx
