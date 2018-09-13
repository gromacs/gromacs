/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Tests for gmx nmr.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/gmxana/gmx_ana.h"

#include "testutils/cmdlinetest.h"
#include "testutils/stdiohelper.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

namespace
{

using gmx::test::XvgMatch;

class NmrTest : public gmx::test::CommandLineTestBase
{
    public:
        void runTest()
        {
            auto &cmdline = commandLine();

            setInputFile("-s", "nmr.tpr");
            setInputFile("-f", "nmr.edr");

            gmx::test::StdioTestHelper stdioHelper(&fileManager());
            stdioHelper.redirectStringToStdin("-1\n");
            ASSERT_EQ(0, gmx_nmr(cmdline.argc(), cmdline.argv()));

            checkOutputFiles();
        }
};

TEST_F(NmrTest, WritesOrientA)
{
    setOutputFile("-ora", "orienta.xvg", XvgMatch());
    runTest();
}

TEST_F(NmrTest, WritesOrientT)
{
    setOutputFile("-ort", "orientt.xvg", XvgMatch());
    runTest();
}

TEST_F(NmrTest, WritesOrientDeviationA)
{
    setOutputFile("-oda", "orideva.xvg", XvgMatch());
    runTest();
}

TEST_F(NmrTest, WritesOrientDeviationT)
{
    setOutputFile("-odt", "oriedevt.xvg", XvgMatch());
    runTest();
}

TEST_F(NmrTest, WritesOrientDeviationR)
{
    setOutputFile("-odr", "oriedevr.xvg", XvgMatch());
    runTest();
}



} // namespace
