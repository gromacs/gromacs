/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Tests for functionality of the "trajectory" trajectory analysis module.
 *
 * Line coverage for the module code is good, but due to missing input data
 * with velocities or forces, the output of these quantities is not really
 * tested.  But the values are only computed in the selection engine, and
 * at least the low-level computation is tested there.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/trajectory.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"
#include "testutils/textblockmatchers.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::NoTextMatch;

/********************************************************************
 * Tests for gmx::analysismodules::Trajectory.
 */

//! Test fixture for the select analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::TrajectoryInfo> TrajectoryModuleTest;

TEST_F(TrajectoryModuleTest, BasicTest)
{
    const char* const cmdline[] = { "trajectory", "-select", "resnr 1", "resnr 3" };
    setTopology("simple.gro");
    setTrajectory("simple.gro");
    setOutputFile("-ox", "coord.xvg", NoTextMatch());
    includeDataset("x");
    runTest(CommandLine(cmdline));
}

TEST_F(TrajectoryModuleTest, PlotsXOnly)
{
    const char* const cmdline[] = { "trajectory", "-select", "resnr 1", "resnr 3", "-x" };
    setTopology("simple.gro");
    setTrajectory("simple.gro");
    setOutputFile("-ox", "coord.xvg", NoTextMatch());
    includeDataset("x");
    runTest(CommandLine(cmdline));
}

TEST_F(TrajectoryModuleTest, HandlesNoVelocities)
{
    const char* const cmdline[] = {
        "trajectory",
        "-select",
        "resnr 1",
        "resnr 3",
    };
    setTopology("simple.gro");
    setTrajectory("simple.gro");
    setOutputFile("-ov", "vel.xvg", NoTextMatch());
    includeDataset("v");
    runTest(CommandLine(cmdline));
}

TEST_F(TrajectoryModuleTest, HandlesNoForces)
{
    const char* const cmdline[] = {
        "trajectory",
        "-select",
        "resnr 1",
        "resnr 3",
    };
    setTopology("simple.gro");
    setTrajectory("simple.gro");
    setOutputFile("-of", "force.xvg", NoTextMatch());
    includeDataset("f");
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace test
} // namespace gmx
