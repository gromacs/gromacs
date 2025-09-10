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
 * \brief Tests for gmx saltbr.
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

class GmxSaltbrTest : public GmxAnaTestBase
{
    int gmxTool(int argc, char* argv[]) const override { return gmx_saltbr(argc, argv); }
};

TEST_F(GmxSaltbrTest, BasicOutputWorks)
{
    // saltbr creates huge output files, so we check a tiny trajectory of 2 molecules
    const std::string fileNameBase = "spc2-traj";
    TprAndFileManager tprFileHandle(fileNameBase);

    commandLine().addOption("-f", TestFileManager::getInputFilePath(fileNameBase + ".xtc"));
    commandLine().addOption("-s", tprFileHandle.tprName());

    const FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(1.0, 1e-4);
    setOutputFileWithGeneratedName("min-min.xvg", "min-min.xvg", XvgMatch().tolerance(tolerance));
    setOutputFileWithGeneratedName("plus-min.xvg", "plus-min.xvg", XvgMatch().tolerance(tolerance));
    setOutputFileWithGeneratedName("plus-plus.xvg", "plus-plus.xvg", XvgMatch().tolerance(tolerance));

    selectGroups({ "System" });
    runAndCheckResults();
}

} // namespace
} // namespace test
} // namespace gmx
