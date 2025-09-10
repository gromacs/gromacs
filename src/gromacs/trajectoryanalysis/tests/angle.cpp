/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * Tests for functionality of the "angle" trajectory analysis module.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/angle.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;

/********************************************************************
 * Tests for gmx::analysismodules::Angle.
 */

//! Test fixture for the angle analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::AngleInfo> AngleModuleTest;

TEST_F(AngleModuleTest, ComputesSimpleAngles)
{
    const char* const cmdline[] = {
        "angle", "-g1", "angle", "-group1", "resname RA1 RA2 and name A1 A2 A3", "-binw", "60"
    };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, ComputesDihedrals)
{
    const char* const cmdline[] = {
        "angle", "-g1", "dihedral", "-group1", "resname RD1 RD2 RD3 and name A1 A2 A3 A4",
        "-binw", "120"
    };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorPairAngles)
{
    const char* const cmdline[] = { "angle",
                                    "-g1",
                                    "vector",
                                    "-group1",
                                    "resname RV1 RV2 and name A1 A2",
                                    "-g2",
                                    "vector",
                                    "-group2",
                                    "resname RV3 RV4 and name A1 A2",
                                    "-binw",
                                    "60" };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorPlanePairAngles)
{
    const char* const cmdline[] = { "angle",
                                    "-g1",
                                    "vector",
                                    "-group1",
                                    "resname RV1 RV2 and name A1 A2",
                                    "-g2",
                                    "plane",
                                    "-group2",
                                    "resname RP1 RP2 and name A1 A2 A3",
                                    "-binw",
                                    "60" };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, ComputesPlaneZAxisAngles)
{
    const char* const cmdline[] = {
        "angle", "-g1", "plane", "-group1", "resname RP1 RP2 and name A1 A2 A3",
        "-g2",   "z",   "-binw", "60"
    };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorSphereNormalZAxisAngles)
{
    const char* const cmdline[] = {
        "angle", "-g1",     "vector",  "-group1",           "resname RV1 RV2 and name A1 A2",
        "-g2",   "sphnorm", "-group2", "cog of resname RS", "-binw",
        "60"
    };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorTimeZeroAngles)
{
    const char* const cmdline[] = {
        "angle", "-g1", "vector", "-group1", "resname RV1 RV2 RV3 RV4 and name A1 A2",
        "-g2",   "t0",  "-binw",  "60"
    };
    setTopology("angle.gro");
    setTrajectory("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, ComputesMultipleAngles)
{
    const char* const cmdline[] = { "angle",
                                    "-g1",
                                    "vector",
                                    "-group1",
                                    "resname RV1 RV2 and name A1 A2",
                                    "resname RV3 RV4 and name A1 A2",
                                    "-g2",
                                    "plane",
                                    "-group2",
                                    "resname RP1 RP2 and name A1 A2 A3",
                                    "resname RP1 RP2 and name A1 A2 A3",
                                    "-binw",
                                    "60" };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, HandlesDynamicSelections)
{
    const char* const cmdline[] = {
        "angle", "-g1", "angle", "-group1", "resname RA1 RA2 and name A1 A2 A3 and z < 0.5",
        "-binw", "60"
    };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, HandlesOneVsMultipleVectorAngles)
{
    const char* const cmdline[] = { "angle",
                                    "-g1",
                                    "vector",
                                    "-group1",
                                    "resname RV1 RV2 and name A1 A2",
                                    "-g2",
                                    "vector",
                                    "-group2",
                                    "resname RV3 and name A1 A2",
                                    "-binw",
                                    "60" };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(AngleModuleTest, HandlesOneVsMultipleVectorGroupsAngles)
{
    const char* const cmdline[] = { "angle",
                                    "-g1",
                                    "vector",
                                    "-group1",
                                    "resname RV2 and name A1 A2",
                                    "resname RV3 RV4 and name A1 A2",
                                    "-g2",
                                    "plane",
                                    "-group2",
                                    "resname RP1 RP2 and name A1 A2 A3",
                                    "-binw",
                                    "60" };
    setTopology("angle.gro");
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace test
} // namespace gmx
