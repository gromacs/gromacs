/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2016,2017,2018, by the GROMACS development team, led by
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
 * Tests for gmx trjconv.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/gmxana/gmx_ana.h"

#include "testutils/cmdlinetest.h"
#include "testutils/stdiohelper.h"
#include "testutils/textblockmatchers.h"

namespace
{

class TrjconvWithIndexGroupSubset : public gmx::test::CommandLineTestBase,
                                    public ::testing::WithParamInterface<const char *>
{
    public:
        void runTest(const char *fileName)
        {
            auto &cmdline = commandLine();

            setInputFile("-s", "spc2.gro");
            setInputFile("-f", fileName);
            setInputFile("-n", "spc2.ndx");
            setOutputFile("-o", "spc-traj.tng", gmx::test::NoTextMatch());

            gmx::test::StdioTestHelper stdioHelper(&fileManager());
            stdioHelper.redirectStringToStdin("SecondWaterMolecule\n");

            /* TODO Ideally, we would then check that the output file
               has only 3 of the 6 atoms (which it does), but the
               infrastructure for doing that automatically is still
               being built. This would also fix the TODO below. */
            ASSERT_EQ(0, gmx_trjconv(cmdline.argc(), cmdline.argv()));
        }
};
/* TODO These tests are actually not very effective, because trjconv
 * can only return 0 or exit via gmx_fatal() (which currently also
 * exits the test binary). So, "no XDR/TNG support in the binary"
 * leads to the reading test appearing to pass. Still, no fatal error
 * and no segfault is useful information while modifying such code. */

TEST_P(TrjconvWithIndexGroupSubset, WithDifferentInputFormats)
{
    runTest(GetParam());
}

/*! \brief Helper array of input files present in the source repo
 * database. These all have two identical frames of two SPC water
 * molecules, which were generated via trjconv from the .gro
 * version. */
const char *const trajectoryFileNames[] = {
    "spc2-traj.trr",
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.xtc",
    "spc2-traj.gro",
    "spc2-traj.pdb",
    "spc2-traj.g96"
};

INSTANTIATE_TEST_CASE_P(NoFatalErrorWhenWritingFrom,
                        TrjconvWithIndexGroupSubset,
                            ::testing::ValuesIn(trajectoryFileNames));

} // namespace
