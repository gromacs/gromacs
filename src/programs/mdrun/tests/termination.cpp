/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Tests for the mdrun termination functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <gtest/gtest.h>

#include "testutils/testfilemanager.h"

#include "moduletest.h"
#include "terminationhelper.h"

namespace gmx
{
namespace test
{

//! Build a simple .mdp file
void organizeMdpFile(SimulationRunner *runner)
{
    // Make sure -maxh has a chance to propagate
    runner->useStringAsMdpFile("nsteps = 100\n"
                               "tcoupl = v-rescale\n"
                               "tc-grps = System\n"
                               "tau-t = 1\n"
                               "ref-t = 298\n");
}

//! Convenience typedef
typedef MdrunTestFixture MdrunTerminationTest;

TEST_F(MdrunTerminationTest, WritesCheckpointAfterMaxhTerminationAndThenRestarts)
{
    CommandLine       mdrunCaller;
    mdrunCaller.append("mdrun");
    TerminationHelper helper(&fileManager_, &mdrunCaller, &runner_);

    organizeMdpFile(&runner_);
    EXPECT_EQ(0, runner_.callGrompp());

    helper.runFirstMdrun(runner_.cptFileName_);
    helper.runSecondMdrun();
}

} // namespace
} // namespace
