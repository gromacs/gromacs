/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Tests for functionality of the "pairdist" trajectory analysis module.
 *
 * These tests test the basic functionality of the tool itself, but currently
 * the following are missing:
 *  - Tests related to -odg output.  This would require a full tpr file, and
 *    some investigation on what kind of tpr it should be to produce reasonable
 *    output.
 *  - Tests for the X axes in the area per atom/residue plots.  These could be
 *    added once better X axes are implemented.
 *  - Tests for XVG labels.  This is a limitation of the current testing
 *    framework.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/pairdist.h"

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;

/********************************************************************
 * Tests for gmx::analysismodules::PairDistance.
 */

//! Test fixture for the select analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::PairDistanceInfo>
    PairDistanceModuleTest;

TEST_F(PairDistanceModuleTest, ComputesAllDistances)
{
    const char *const cmdline[] = {
        "pairdist",
        "-ref", "resindex 1", "-refgrouping", "none",
        "-sel", "resindex 3", "-selgrouping", "none"
    };
    setTopology("simple.gro");
    setOutputFileNoTest("-o", "xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesAllDistancesWithCutoff)
{
    const char *const cmdline[] = {
        "pairdist",
        "-ref", "resindex 1", "-refgrouping", "none",
        "-sel", "resindex 3", "-selgrouping", "none",
        "-cutoff", "1.5"
    };
    setTopology("simple.gro");
    setOutputFileNoTest("-o", "xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesMinDistanceWithCutoff)
{
    const char *const cmdline[] = {
        "pairdist",
        "-ref", "resindex 1",
        "-sel", "resindex 3",
        "-cutoff", "1.5"
    };
    setTopology("simple.gro");
    setOutputFileNoTest("-o", "xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesMaxDistance)
{
    const char *const cmdline[] = {
        "pairdist",
        "-ref", "resindex 1",
        "-sel", "resindex 3",
        "-type", "max"
    };
    setTopology("simple.gro");
    setOutputFileNoTest("-o", "xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesMaxDistanceWithCutoff)
{
    const char *const cmdline[] = {
        "pairdist",
        "-ref", "resindex 1",
        "-sel", "resindex 3",
        "-cutoff", "1.5", "-type", "max"
    };
    setTopology("simple.gro");
    setOutputFileNoTest("-o", "xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesGroupedMinDistanceWithCutoff)
{
    const char *const cmdline[] = {
        "pairdist",
        "-ref", "resindex 1 to 2", "-refgrouping", "res",
        "-sel", "resindex 3 to 5", "-selgrouping", "res",
        "-cutoff", "2.5"
    };
    setTopology("simple.gro");
    setOutputFileNoTest("-o", "xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesGroupedMaxDistanceWithCutoff)
{
    const char *const cmdline[] = {
        "pairdist",
        "-ref", "resindex 1 to 2", "-refgrouping", "res",
        "-sel", "resindex 3 to 5", "-selgrouping", "res",
        "-cutoff", "3.5", "-type", "max"
    };
    setTopology("simple.gro");
    setOutputFileNoTest("-o", "xvg");
    runTest(CommandLine(cmdline));
}

} // namespace
