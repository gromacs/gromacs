/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * This implements basic initial constrains test (using single-rank mdrun)
 *
 * This test checks that the coordinates from file, which only satisfy
 * the constraints up to gro precision, are constrained correctly and that
 * the initial velocity of the center of mass does not contribute to the
 * kinetic energy at step 0..
 * It runs the input system for 1 step (no continuation), and compares the total energy.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! This type holds input integrators. Now it's holding names, but ei* enum values from md_enums.h could be used instead.
using EnergyIntegratorType = const char*;

//! Test fixture parametrized on integrators
class InitialConstraintsTest :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<EnergyIntegratorType>
{
};

TEST_P(InitialConstraintsTest, Works)
{
    const int         nsteps     = 1;
    const float       timestep   = 0.001;
    const auto*       integrator = GetParam();
    const std::string integratorName(integrator);
    SCOPED_TRACE("Integrating with " + integratorName);
    const std::string theMdpFile = formatString(
            "nstcalcenergy           = 1\n"
            "nstenergy               = 1\n"
            "comm-mode               = linear\n"
            "continuation            = no\n"
            "constraints             = h-bonds\n"
            "lincs_iter              = 2\n"
            "verlet-buffer-tolerance = 1e-4\n"
            "nsttcouple              = 1\n" // for md-vv-avek
            "nstpcouple              = 1\n" // for md-vv-avek
            "integrator              = %s\n"
            "nsteps                  = %d\n"
            "dt                      = %f\n",
            integratorName.c_str(),
            nsteps,
            timestep);

    runner_.useStringAsMdpFile(theMdpFile);

    const std::string inputFile = "spc-and-methanol";
    runner_.useTopGroAndNdxFromDatabase(inputFile);
    EXPECT_EQ(0, runner_.callGrompp());

    runner_.edrFileName_ = fileManager_.getTemporaryFilePath(inputFile + ".edr").string();
    ASSERT_EQ(0, runner_.callMdrun());

    auto energyReader =
            openEnergyFileToReadTerms(runner_.edrFileName_, { "Total Energy", "Kinetic En." });
    real totalEnergy = 0.0, prevTotalEnergy = 0.0;
    auto tolerance = ulpTolerance(0); // The real value is set below from starting kinetic energy
    for (int i = 0; i <= nsteps; i++)
    {
        EnergyFrame frame = energyReader->frame();
        prevTotalEnergy   = totalEnergy;
        totalEnergy       = frame.at("Total Energy");
        if (i == 0)
        {
            // We set the tolerance for total energy based on magnitude of kinetic energy.
            // The reason is that the other total energy component, the potential energy, can in theory have whatever magnitude.
            const real startingKineticEnergy = frame.at("Kinetic En.");
            tolerance = relativeToleranceAsFloatingPoint(startingKineticEnergy, 1e-5);
        }
        else
        {
            EXPECT_REAL_EQ_TOL(totalEnergy, prevTotalEnergy, tolerance);
        }
    }
}

//! Integrators with energy conservation to test
const EnergyIntegratorType c_integratorsToTest[] = { "md", "md-vv", "md-vv-avek" };

INSTANTIATE_TEST_SUITE_P(Checking, InitialConstraintsTest, ::testing::ValuesIn(c_integratorsToTest));

} // namespace
} // namespace test
} // namespace gmx
