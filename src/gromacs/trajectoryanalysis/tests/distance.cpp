/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Tests for functionality of the "distance" trajectory analysis module.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/distance.h"

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;

/********************************************************************
 * Tests for gmx::analysismodules::Distance.
 */

//! Test fixture for the angle analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::DistanceInfo>
    DistanceModuleTest;

TEST_F(DistanceModuleTest, ComputesDistances)
{
    const char *const cmdline[] = {
        "distance",
        "-select", "atomname S1 S2",
        "-len", "2", "-binw", "0.5"
    };
    setTopology("simple.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(DistanceModuleTest, ComputesMultipleDistances)
{
    const char *const cmdline[] = {
        "distance",
        "-select", "atomname S1 S2",
        "resindex 1 to 4 and atomname CB merge resindex 2 to 5 and atomname CB",
        "-len", "2", "-binw", "0.5"
    };
    setTopology("simple.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(DistanceModuleTest, HandlesDynamicSelections)
{
    const char *const cmdline[] = {
        "distance",
        "-select", "atomname S1 S2 and res_cog x < 2.8",
        "-len", "2", "-binw", "0.5"
    };
    setTopology("simple.gro");
    runTest(CommandLine(cmdline));
}

} // namespace
