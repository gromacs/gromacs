/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests for gmx traj.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxana/gmx_ana.h"

#include "testutils/cmdlinetest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/stdiohelper.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

namespace gmx
{
namespace test
{
namespace
{

class GmxTraj : public gmx::test::CommandLineTestBase, public ::testing::WithParamInterface<const char*>
{
public:
    void runTest(const char* fileName)
    {
        auto& cmdline = commandLine();
        setInputFile("-s", "spc2.gro");
        setInputFile("-f", fileName);
        setOutputFile("-ox", "spc2.xvg", gmx::test::NoTextMatch());

        gmx::test::StdioTestHelper stdioHelper(&fileManager());
        stdioHelper.redirectStringToStdin("0\n");

        ASSERT_EQ(0, gmx_traj(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }
};

/* TODO These tests are actually not very effective, because gmx-traj
 * can only return 0 or exit via gmx_fatal() (which currently also
 * exits the test binary). So, "no XDR/TNG support in the binary"
 * leads to the reading test appearing to pass. Still, no fatal error
 * and no segfault is useful information while modifying such code. */

TEST_P(GmxTraj, WithDifferentInputFormats)
{
    runTest(GetParam());
}

TEST_P(GmxTraj, RotationalKineticEnergy)
{
    setOutputFile("-ekr",
                  "spc2_ekr.xvg",
                  gmx::test::XvgMatch().tolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-4)));
    runTest(GetParam());
}

/*! \brief Helper array of input files present in the source repo
 * database. These all have two identical frames of two SPC water
 * molecules, which were generated via trjconv from the .gro
 * version. */
const char* const trajectoryFileNames[] = { "spc2-traj.trr", "spc2-traj.xtc", "spc2-traj.gro",
                                            "spc2-traj.pdb", "spc2-traj.g96",
#if GMX_USE_TNG
                                            "spc2-traj.tng"
#endif
};

INSTANTIATE_TEST_SUITE_P(NoFatalErrorWhenWritingFrom, GmxTraj, ::testing::ValuesIn(trajectoryFileNames));

} // namespace
} // namespace test
} // namespace gmx
