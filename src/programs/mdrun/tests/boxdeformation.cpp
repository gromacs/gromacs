/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Test for continous box deformation
 *
 * Testing that one actually obtains the correct stress for a given
 * box deformation rate requires an extreme amount of sampling.
 * Here we mainly check that the system is stable and that the velocity
 * is corrected for atoms that move over PBC.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"

#include "energycomparison.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Test fixture base for simple mdrun systems
 *
 * This test ensures mdrun can run a simulation, reaching reproducible energies.
 * This test will fail when velocities are not corrected when moving atoms
 * by a periodic box vector.
 */
class BoxDeformationTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string>>
{
};

// Check that a system with Ekin=0 reports Ekin=0 under deformation
TEST_F(BoxDeformationTest, flowDoesNotAffectEkin)
{
    // We use an LJ system with cut-off < sigma to have zero forces
    std::string theMdpFile =
            "coulombtype      = reaction-field\n"
            "nstenergy        = 10\n"
            "nstlist          = 10\n"
            "verlet-buffer-tolerance = -1\n"
            "rlist            = 0.33\n"
            "rcoulomb         = 0.3\n"
            "rvdw             = 0.3\n"
            "dt               = 0.002\n"
            "nsteps           = 0\n"
            "gen-vel          = yes\n"
            "gen-temp         = 0\n"
            "deform           = 0 0 0 0 6e-2 0\n"
            "deform-init-flow = yes\n";

    const auto& simulationName = "argon12";

    // Prepare the .tpr file
    {
        CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(theMdpFile);
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Do mdrun
    {
        CommandLine mdrunCaller;
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
        auto relativeTolerance = relativeToleranceAsFloatingPoint(0.01, GMX_DOUBLE ? 1e-6 : 1e-4);
        EnergyTermsToCompare energyTermsToCompare{ { interaction_function[F_EKIN].longname,
                                                     relativeTolerance } };
        TestReferenceData    refData;
        auto checker = refData.rootChecker().checkCompound("Simulation", simulationName);
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
    }
}

// Check that a deformed water box has correct Epot and Ekin
TEST_F(BoxDeformationTest, EnergiesWithinTolerances)
{
    // We force nstlist=10 to trigger put_atoms_in_box more often
    // which requires correcting the velocities with the box deformation.
    // We use a (high) shear rate of 10 (ns^-1).
    std::string theMdpFile =
            "coulombtype      = PME\n"
            "nstenergy        = 10\n"
            "nstlist          = 10\n"
            "verlet-buffer-tolerance = -1\n"
            "rlist            = 0.93\n"
            "rcoulomb         = 0.9\n"
            "rvdw             = 0.9\n"
            "pme-order        = 4\n"
            "fourier-spacing  = 0.12\n"
            "dt               = 0.002\n"
            "nsteps           = 20\n"
            "deform           = 0 0 0 0 1.86e-2 0\n"
            "deform-init-flow = yes\n";

    const auto& simulationName = "spc216";

    // Prepare the .tpr file
    {
        CommandLine caller;
        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(theMdpFile);
        EXPECT_EQ(0, runner_.callGrompp(caller));
    }
    // Do mdrun
    {
        CommandLine mdrunCaller;
        ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));
        auto relativeTolerance = relativeToleranceAsFloatingPoint(1, GMX_DOUBLE ? 1e-6 : 1e-4);
        EnergyTermsToCompare energyTermsToCompare{
            { { interaction_function[F_EPOT].longname, relativeTolerance },
              { interaction_function[F_EKIN].longname, relativeTolerance } }
        };
        TestReferenceData refData;
        auto checker = refData.rootChecker().checkCompound("Simulation", simulationName);
        checkEnergiesAgainstReferenceData(runner_.edrFileName_, energyTermsToCompare, &checker);
    }
}

} // namespace
} // namespace test
} // namespace gmx
