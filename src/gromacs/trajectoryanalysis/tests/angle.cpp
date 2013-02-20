/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Tests for functionality of the "angle" trajectory analysis module.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include <gtest/gtest.h>

#include "gromacs/trajectoryanalysis/modules/angle.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;

/********************************************************************
 * Tests for gmx::analysismodules::Angle.
 */

//! Test fixture for the angle analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::Angle>
    AngleModuleTest;

TEST_F(AngleModuleTest, ComputesSimpleAngles)
{
    const char *const cmdline[] = {
        "angle",
        "-g1", "angle", "-group1", "resname RA1 RA2 and name A1 A2 A3"
    };
    setTopology("angle.gro");
    runTest(CommandLine::create(cmdline));
}

TEST_F(AngleModuleTest, ComputesDihedrals)
{
    const char *const cmdline[] = {
        "angle",
        "-g1", "dihedral", "-group1", "resname RD1 RD2 RD3 and name A1 A2 A3 A4"
    };
    setTopology("angle.gro");
    runTest(CommandLine::create(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorPairAngles)
{
    const char *const cmdline[] = {
        "angle",
        "-g1", "vector", "-group1", "resname RV1 RV2 and name A1 A2",
        "-g2", "vector", "-group2", "resname RV3 RV4 and name A1 A2"
    };
    setTopology("angle.gro");
    runTest(CommandLine::create(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorPlanePairAngles)
{
    const char *const cmdline[] = {
        "angle",
        "-g1", "vector", "-group1", "resname RV1 RV2 and name A1 A2",
        "-g2", "plane",  "-group2", "resname RP1 RP2 and name A1 A2 A3"
    };
    setTopology("angle.gro");
    runTest(CommandLine::create(cmdline));
}

TEST_F(AngleModuleTest, ComputesPlaneZAxisAngles)
{
    const char *const cmdline[] = {
        "angle",
        "-g1", "plane", "-group1", "resname RP1 RP2 and name A1 A2 A3",
        "-g2", "z"
    };
    setTopology("angle.gro");
    runTest(CommandLine::create(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorSphereNormalZAxisAngles)
{
    const char *const cmdline[] = {
        "angle",
        "-g1", "vector",  "-group1", "resname RV1 RV2 and name A1 A2",
        "-g2", "sphnorm", "-group2", "cog of resname RS"
    };
    setTopology("angle.gro");
    runTest(CommandLine::create(cmdline));
}

TEST_F(AngleModuleTest, ComputesVectorTimeZeroAngles)
{
    const char *const cmdline[] = {
        "angle",
        "-g1", "vector", "-group1", "resname RV1 RV2 RV3 RV4 and name A1 A2",
        "-g2", "t0"
    };
    setTopology("angle.gro");
    setTrajectory("angle.gro");
    runTest(CommandLine::create(cmdline));
}

} // namespace
