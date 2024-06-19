/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Simple tests for the mdrun functionality.
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
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/trajectoryreader.h"
#include "testutils/xvgtest.h"

#include "energycomparison.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Database of enerngy tolerances for MD integrator on the various systems. */
std::unordered_map<std::string, FloatingPointTolerance> energyToleranceForSystem_g = {
    { { "angles1", relativeToleranceAsFloatingPoint(1, 1e-4) } }
};

/*! \brief Database of pressure
   tolerances for MD integrator on the various systems. */
std::unordered_map<std::string, FloatingPointTolerance> pressureToleranceForSystem_g = {
    { { "angles1", relativeToleranceAsFloatingPoint(1, 1e-4) } }
};

//! Helper type
using MdpField = MdpFieldValues::value_type;

/*! \brief Test fixture base for simple mdrun systems
 *
 * This test ensures mdrun can run a simulation, reaching
 * reproducible energies.
 *
 * The choices for tolerance are arbitrary but sufficient. */
class SimpleMdrunTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string>>
{
};

TEST_P(SimpleMdrunTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto integrator     = std::get<1>(params);
    SCOPED_TRACE(formatString("Comparing simple mdrun for '%s'", simulationName.c_str()));

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
    mdpFieldValues["nsteps"]        = "50";
    mdpFieldValues["nstfout"]       = "4";
    mdpFieldValues["constraints"]   = "none";
    mdpFieldValues["nstcalcenergy"] = "4";
    mdpFieldValues["coulombtype"]   = "Cut-off";
    mdpFieldValues["vdwtype"]       = "Cut-off";

    // Prepare the .tpr file
    {
        CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Do mdrun
    {
        CommandLine mdrunCaller;
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
        EnergyTermsToCompare energyTermsToCompare{ {
                { interaction_function[F_EPOT].longname, energyToleranceForSystem_g.at(simulationName) },
                { interaction_function[F_EKIN].longname, energyToleranceForSystem_g.at(simulationName) },
                { interaction_function[F_PRES].longname, pressureToleranceForSystem_g.at(simulationName) },
        } };
        TestReferenceData    refData;
        auto                 checker = refData.rootChecker()
                               .checkCompound("Simulation", simulationName)
                               .checkCompound("Mdrun", integrator);
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
        // Now check the forces
        TrajectoryFrameReader reader(runner_.fullPrecisionTrajectoryFileName_);
        checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(1, 1e-4));
        do
        {
            auto frame = reader.frame();
            auto force = frame.f();
            int  atom  = 0;
            for (const auto& f : force)
            {
                std::string forceName = frame.frameName() + " F[" + toString(atom) + "]";

                checker.checkVector(f, forceName.c_str());
                atom++;
            }
        } while (reader.readNextFrame());
    }
}

// The time for OpenCL kernel compilation means these tests might time
// out. If that proves to be a problem, these can be disabled for
// OpenCL builds. However, once that compilation is cached for the
// lifetime of the whole test binary process, these tests should run in
// such configurations.
#if GMX_DOUBLE

//! Containers of systems to test.
//! \{
std::vector<std::string> systemsToTest_g = { "angles1" };
std::vector<std::string> md_g            = { "md", "md-vv" };
//! \}

INSTANTIATE_TEST_SUITE_P(Angles1,
                         SimpleMdrunTest,
                         ::testing::Combine(::testing::ValuesIn(systemsToTest_g),
                                            ::testing::ValuesIn(md_g)));
#else
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(SimpleMdrunTest);
#endif
} // namespace
} // namespace test
} // namespace gmx
