/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
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
 * Tests for gmx::analysismodules::PairDistance.
 */

//! Test fixture for the select analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::PairDistanceInfo> PairDistanceModuleTest;

TEST_F(PairDistanceModuleTest, ComputesAllDistances)
{
    const char* const cmdline[] = { "pairdist",     "-ref",         "resindex 1",
                                    "-refgrouping", "none",         "-sel",
                                    "resindex 3",   "-selgrouping", "none" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesAllDistancesWithCutoff)
{
    const char* const cmdline[] = { "pairdist", "-ref",    "resindex 1", "-refgrouping",
                                    "none",     "-sel",    "resindex 3", "-selgrouping",
                                    "none",     "-cutoff", "1.5" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesMinDistanceWithCutoff)
{
    const char* const cmdline[] = { "pairdist",   "-ref",    "resindex 1", "-sel",
                                    "resindex 3", "-cutoff", "1.5" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesMaxDistance)
{
    const char* const cmdline[] = { "pairdist",   "-ref",  "resindex 1", "-sel",
                                    "resindex 3", "-type", "max" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesMaxDistanceWithCutoff)
{
    const char* const cmdline[] = { "pairdist", "-ref", "resindex 1", "-sel", "resindex 3",
                                    "-cutoff",  "1.5",  "-type",      "max" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesGroupedMinDistanceWithCutoff)
{
    const char* const cmdline[] = { "pairdist",        "-ref",         "resindex 1 to 2",
                                    "-refgrouping",    "res",          "-sel",
                                    "resindex 3 to 5", "-selgrouping", "res",
                                    "-cutoff",         "2.5" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, ComputesGroupedMaxDistanceWithCutoff)
{
    const char* const cmdline[] = { "pairdist",
                                    "-ref",
                                    "resindex 1 to 2",
                                    "-refgrouping",
                                    "res",
                                    "-sel",
                                    "resindex 3 to 5",
                                    "-selgrouping",
                                    "res",
                                    "-cutoff",
                                    "3.5",
                                    "-type",
                                    "max" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, CoordinateSelectionIsNotOverwritten)
{
    const char* const cmdline[] = { "pairdist", "-ref", "[0.0, 1.5, 2.9]", "-sel", "resindex 3",
                                    "-type",    "max" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(PairDistanceModuleTest, CoordinateSelectionIsNotOverwrittenWithExplicitGroup)
{
    const char* const cmdline[] = { "pairdist", "-ref", "[0.0, 1.5, 2.9]", "-sel", "resindex 3",
                                    "-type",    "max",  "-refgrouping",    "res" };
    setTopology("simple.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace test
} // namespace gmx
