/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests for basic grompp functionality
 *
 * \todo Refactor SimulationRunner to split off SimulationPreparer, so
 * that integration tests of grompp can stand apart from tests of
 * mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"
#include "testutils/testexceptions.h"

#include "moduletest.h"

namespace
{

//! Test fixture for grompp
class GromppTest : public gmx::test::MdrunTestFixture
{
public:
    //! Execute the trajectory writing test
    void runTest()
    {
        runner_.useTopGroAndNdxFromDatabase("spc-and-methanol");
        EXPECT_EQ(0, runner_.callGrompp());
    }
};

/* This test ensures that an empty .mdp file (ie. all default values) works. */
TEST_F(GromppTest, EmptyMdpFileWorks)
{
    runner_.useEmptyMdpFile();
    runTest();
}

/* Test for making sure grompp can handle simulated annealing data */
TEST_F(GromppTest, SimulatedAnnealingWorks)
{
    runner_.useStringAsMdpFile(
            "annealing = periodic\n"
            "annealing-npoints = 4\n"
            "annealing-time = 0 2 4 6\n"
            "annealing-temp = 298 320 320 298\n");
    runTest();
}

TEST_F(GromppTest, SimulatedAnnealingWorksWithMultipleGroups)
{
    runner_.useStringAsMdpFile(
            "tc-grps = Methanol SOL\n"
            "tau-t = 0.1 0.1\n"
            "ref_t = 298 298\n"
            "annealing = single periodic\n"
            "annealing-npoints = 3 4\n"
            "annealing-time = 0 3 6 0 2 4 6\n"
            "annealing-temp = 298 280 270 298 320 320 298\n");
    runTest();
}

/* Test for making sure that we need to set maxwarn correctly for grompp to run */
TEST_F(GromppTest, DeathTestHandlesNoMaxwarnError)
{
    runner_.useStringAsMdpFile(
            "Tcoupl  = Berendsen\n"
            "tc-grps = System\n"
            "tau-t   = 0.1\n"
            "ref_t   = 298\n");
    GMX_EXPECT_DEATH_IF_SUPPORTED(runTest(), "");
}

TEST_F(GromppTest, HandlesMaxwarn)
{
    runner_.useStringAsMdpFile(
            "Tcoupl  = Berendsen\n"
            "tc-grps = System\n"
            "tau-t   = 0.1\n"
            "ref_t   = 298\n");
    runner_.setMaxWarn(1);
    runTest();
}

TEST_F(GromppTest, MaxwarnShouldBePositive)
{
    runner_.useStringAsMdpFile(
            "Tcoupl  = Berendsen\n"
            "tc-grps = System\n"
            "tau-t   = 0.1\n"
            "ref_t   = 298\n");
    runner_.setMaxWarn(-1);
    EXPECT_THROW_GMX(runTest(), gmx::InconsistentInputError);
}

#if HAVE_MUPARSER

TEST_F(GromppTest, ValidTransformationCoord)
{
    const char* inputMdpFile[] = {
        "pull = yes",
        "pull-ncoords = 2",
        "pull-ngroups = 2",
        "pull-group1-name = SOL",
        "pull-group2-name = Methanol",
        "pull-coord1-geometry = distance",
        "pull-coord1-groups = 1 2",
        "pull-coord2-geometry = transformation",
        "pull-coord2-expression = x1", // Valid expression
    };
    runner_.useStringAsMdpFile(gmx::joinStrings(inputMdpFile, "\n"));
    runTest();
}

TEST_F(GromppTest, InvalidTransformationCoord)
{
    const char* inputMdpFile[] = {
        "pull = yes",
        "pull-ncoords = 2",
        "pull-ngroups = 2",
        "pull-group1-name = SOL",
        "pull-group2-name = Methanol",
        "pull-coord1-geometry = distance",
        "pull-coord1-groups = 1 2",
        "pull-coord2-geometry = transformation",
        "pull-coord2-expression = x2", // Invalid expression -> evaluation should fail
    };
    runner_.useStringAsMdpFile(gmx::joinStrings(inputMdpFile, "\n"));
    ASSERT_THROW(runTest(), gmx::InconsistentInputError);
    done_inputrec_strings(); // This allows grompp to be called again in another test
}

TEST_F(GromppTest, RejectCRescaleAndAnisotropic)
{
    const char* inputMdpFile[] = { "integrator              = md",
                                   "nsteps                  = 1",
                                   "tcoupl                  = v-rescale",
                                   "tc-grps                 = System ",
                                   "ref-t                   = 300  ",
                                   "tau-t                   = 0.1 ",
                                   "pcoupl                  = C-rescale",
                                   "pcoupltype              = anisotropic",
                                   "compressibility         = 1.0 1.0 1.0 1.0 1.0 1.0",
                                   "ref-p                   = 1.0 1.0 1.0 0.0 0.0 0.0" };
    runner_.useStringAsMdpFile(gmx::joinStrings(inputMdpFile, "\n"));
    GMX_EXPECT_DEATH_IF_SUPPORTED(runTest(), "C-rescale does not support pressure coupling type");
}
#endif // HAVE_MUPARSER


} // namespace
