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
 * Tests for gmx traj
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gromacs/gmxana/gmx_ana.h"
#include "testutils/integrationtests.h"
#include "testutils/cmdlinetest.h"

namespace
{

//! Test fixture for gmx traj
typedef gmx::test::IntegrationTestFixture GmxTrajTest;

class GmxTrajTestNoFatalError : public GmxTrajTest
{
    public:
        GmxTrajTestNoFatalError() : groFileName(fileManager_.getInputFilePath("spc2.gro")),
                                    xvgFileName(fileManager_.getTemporaryFilePath("spc2.xvg"))
        {
        }

        int runTest(const char *fileName)
        {
            gmx::test::CommandLine caller;
            caller.append("traj");

            caller.addOption("-s",  groFileName);
            caller.addOption("-ox", xvgFileName);

            std::string inputTrajectoryFileName = fileManager_.getInputFilePath(fileName).c_str();
            caller.addOption("-f", inputTrajectoryFileName);

            redirectStringToStdin("0\n");

            return gmx_traj(caller.argc(), caller.argv());
        }

        std::string groFileName;
        std::string xvgFileName;
};

/* There are fancy ways to use "value-parameterized tests" to reduce
 * the duplication below, but for the handful of file formats we care
 * about, the code would probably get longer (and need fancy template
 * stuff, too).
 */

/* TODO These tests are actually not very effective, because gmx-traj
 * can only return 0 or exit via gmx_fatal() (which currently also
 * exits the test binary). So, "no XDR/TNG support in the binary"
 * leads to the reading test appearing to pass. Still, no fatal error
 * and no segfault is useful information while modifying such code. */
TEST_F(GmxTrajTestNoFatalError, TRRFile)
{
    runTest("spc2-traj.trr");
}

TEST_F(GmxTrajTestNoFatalError, XTCFile)
{
    runTest("spc2-traj.xtc");
}

TEST_F(GmxTrajTestNoFatalError, GROFile)
{
    runTest("spc2-traj.gro");
}

TEST_F(GmxTrajTestNoFatalError, PDBFile)
{
    runTest("spc2-traj.pdb");
}

TEST_F(GmxTrajTestNoFatalError, G96File)
{
    runTest("spc2-traj.g96");
}

} // namespace
