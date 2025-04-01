/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Tests for gmx sorient.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/tests/gmxanatestbase.h"

#include "testutils/testfilemanager.h"
#include "testutils/tprfilegenerator.h"
#include "testutils/xvgtest.h"

namespace gmx
{
namespace test
{
namespace
{

class GmxSorientTest : public GmxAnaTestBase
{
    int gmxTool(int argc, char* argv[]) const override { return gmx_sorient(argc, argv); }
};

TEST_F(GmxSorientTest, BasicOutputWorks)
{
    const std::string fileNameBase = "alanine_vsite_solvated";
    TprAndFileManager tprFileHandle(fileNameBase);

    commandLine().addOption("-f", TestFileManager::getInputFilePath(fileNameBase + ".xtc"));
    commandLine().addOption("-s", tprFileHandle.tprName());
    commandLine().addOption("-e", 1); // only check 2 frames since this is expensive
    commandLine().addOption("-cbin", 0.5);
    commandLine().addOption("-rbin", 0.25);

    const FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(1.0, 1e-4);
    setOutputFile("-o", "sori.xvg", XvgMatch().tolerance(tolerance));
    setOutputFile("-no", "snor.xvg", XvgMatch().tolerance(tolerance));
    setOutputFile("-ro", "sord.xvg", XvgMatch().tolerance(tolerance));
    setOutputFile("-co", "scum.xvg", XvgMatch().tolerance(tolerance));
    setOutputFile("-rc", "scount.xvg", XvgMatch().tolerance(tolerance));

    selectGroups({ "Protein", "Water" });
    runAndCheckResults();
}

} // namespace
} // namespace test
} // namespace gmx
