/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2019, by the GROMACS development team, led by
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

#include <gtest/gtest.h>

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


} // namespace
