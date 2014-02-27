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
 * Tests for functionality of the "sasa" trajectory analysis module.
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
 * The actual surface area algorithm should be tested separately with a more
 * extensive set of data, but those fit better as a separate set of unit tests.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include <gtest/gtest.h>

#include "gromacs/trajectoryanalysis/modules/sasa.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

#include "moduletest.h"

namespace
{

using gmx::test::CommandLine;

/********************************************************************
 * Tests for gmx::analysismodules::Sas.
 */

//! Test fixture for the select analysis module.
typedef gmx::test::TrajectoryAnalysisModuleTestFixture<gmx::analysismodules::SasaInfo>
    SasaModuleTest;

TEST_F(SasaModuleTest, BasicTest)
{
    const char *const cmdline[] = {
        "sasa",
        "-surface", "all",
        "-output", "name N CA C O H"
    };
    setTopology("lysozyme.gro");
    setOutputFileNoTest("-o", "xvg");
    setOutputFileNoTest("-or", "xvg");
    setOutputFileNoTest("-oa", "xvg");
    setOutputFileNoTest("-tv", "xvg");
    excludeDataset("dgsolv");
    setDatasetTolerance("area", gmx::test::ulpTolerance(8));
    setDatasetTolerance("volume", gmx::test::ulpTolerance(8));
    runTest(CommandLine(cmdline));
}

TEST_F(SasaModuleTest, WritesConnollySurfaceWithSolute)
{
    const char *const cmdline[] = {
        "sasa",
        "-surface", "atomnr 1"
    };
    setTopology("lysozyme.gro");
    setOutputFileNoTest("-o", "xvg");
    setOutputFile("-q", "connolly.pdb");
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
    const char *const cmdline[] = {
        "sasa",
        "-surface", "all",
        "-output", "y > 1.5"
    };
    setTopology("lysozyme.gro");
    setOutputFileNoTest("-o", "xvg");
    setOutputFileNoTest("-or", "xvg");
    setOutputFileNoTest("-oa", "xvg");
    excludeDataset("volume");
    excludeDataset("dgsolv");
    setDatasetTolerance("area", gmx::test::ulpTolerance(8));
    runTest(CommandLine(cmdline));
}

TEST_F(SasaModuleTest, HandlesDynamicCalculationGroup)
{
    const char *const cmdline[] = {
        "sasa",
        "-surface", "y > 1.5",
        "-output", "y > 1.5 and z > 0"
    };
    setTopology("lysozyme.gro");
    setOutputFileNoTest("-o", "xvg");
    setOutputFileNoTest("-or", "xvg");
    setOutputFileNoTest("-oa", "xvg");
    excludeDataset("volume");
    excludeDataset("dgsolv");
    setDatasetTolerance("area", gmx::test::ulpTolerance(8));
    runTest(CommandLine(cmdline));
}

} // namespace
