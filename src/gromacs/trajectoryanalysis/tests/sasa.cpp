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
 * Tests for functionality of the "sasa" trajectory analysis module.
 *
 * These tests test the basic functionality of the tool itself, but currently
 * the tests related to -odg output are missing.  This would require a full tpr
 * file, and some investigation on what kind of tpr it should be to produce
 * reasonable output.
 *
 * The actual surface area algorithm is tested separately in surfacearea.cpp.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/sasa.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::ExactTextMatch;
using gmx::test::NoTextMatch;
using gmx::test::XvgMatch;

/********************************************************************
 * Tests for gmx::analysismodules::Sasa.
 */

//! Test fixture for the `sasa` analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::SasaInfo> SasaModuleTest;

TEST_F(SasaModuleTest, BasicTest)
{
    const char* const cmdline[] = { "sasa", "-surface", "all", "-output", "name N CA C O H" };
    setTopology("lysozyme.gro");
    setOutputFile("-o", ".xvg", XvgMatch().testData(false));
    setOutputFile("-or", ".xvg", XvgMatch());
    setOutputFile("-oa", ".xvg", XvgMatch());
    setOutputFile("-tv", ".xvg", XvgMatch().testData(false));
    excludeDataset("dgsolv");
    setDatasetTolerance("area", gmx::test::ulpTolerance(8));
    setDatasetTolerance("volume", gmx::test::ulpTolerance(8));
    runTest(CommandLine(cmdline));
}

TEST_F(SasaModuleTest, HandlesSelectedResidues)
{
    const char* const cmdline[] = { "sasa", "-surface", "resnr 2 4 to 5 8" };
    setTopology("lysozyme.gro");
    setOutputFile("-o", ".xvg", XvgMatch().testData(false));
    setOutputFile("-or", ".xvg", XvgMatch());
    setOutputFile("-oa", ".xvg", XvgMatch());
    excludeDataset("dgsolv");
    excludeDataset("volume");
    setDatasetTolerance("area", gmx::test::ulpTolerance(8));
    runTest(CommandLine(cmdline));
}

TEST_F(SasaModuleTest, WritesConnollySurfaceWithSolute)
{
    const char* const cmdline[] = { "sasa", "-surface", "atomnr 1" };
    setTopology("lysozyme.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    setOutputFile("-q", "connolly.pdb", ExactTextMatch());
    includeDataset("area");
    runTest(CommandLine(cmdline));
}

// The CONECT records written to the output make this very difficult to test
// without parsing the output file, since by construction, there are a lot of
// dot pairs that are at exactly the same distance from each other in exact
// arithmetic, and depending on rounding etc., different pairs get picked into
// the output for single and double precision.
#if 0
TEST_F(SasaModuleTest, WritesConnollySurfaceWithoutSolute)
{
    const char *const cmdline[] = {
        "sasa",
        "-surface", "atomnr 1",
        "-noprot"
    };
    setTopology("lysozyme.gro");
    setOutputFileNoTest("-o", "xvg");
    setOutputFile("-q", "connolly.pdb");
    includeDataset("area");
    runTest(CommandLine(cmdline));
}
#endif

TEST_F(SasaModuleTest, HandlesDynamicOutputGroup)
{
    const char* const cmdline[] = { "sasa", "-surface", "all", "-output", "y > 1.5" };
    setTopology("lysozyme.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    setOutputFile("-or", ".xvg", NoTextMatch());
    setOutputFile("-oa", ".xvg", NoTextMatch());
    excludeDataset("volume");
    excludeDataset("dgsolv");
    setDatasetTolerance("area", gmx::test::ulpTolerance(8));
    runTest(CommandLine(cmdline));
}

TEST_F(SasaModuleTest, HandlesDynamicCalculationGroup)
{
    const char* const cmdline[] = { "sasa", "-surface", "y > 1.5", "-output", "y > 1.5 and z > 0" };
    setTopology("lysozyme.gro");
    setOutputFile("-o", ".xvg", NoTextMatch());
    setOutputFile("-or", ".xvg", NoTextMatch());
    setOutputFile("-oa", ".xvg", NoTextMatch());
    excludeDataset("volume");
    excludeDataset("dgsolv");
    setDatasetTolerance("area", gmx::test::ulpTolerance(8));
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace test
} // namespace gmx
