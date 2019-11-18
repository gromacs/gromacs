/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Tests for gmx make_ndx.
 *
 * \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#include "gmxpre.h"

#include "gromacs/gmxana/gmx_ana.h"

#include "testutils/cmdlinetest.h"
#include "testutils/stdiohelper.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace
{

class GmxMakeNdx : public gmx::test::CommandLineTestBase
{
public:
    /*! \brief runs the test for the filename prefix \p sysName from the simulation data base
        \param[in]   sysName   file name prefix sysName.g96 must be present in the simulation database
     */
    void runTest(const std::string& sysName)
    {
        auto&       cmdline     = commandLine();
        auto        groFileName = sysName + ".g96";
        std::string ndxFileName = sysName + ".ndx";
        setInputFile("-f", groFileName);
        setOutputFile("-o", ndxFileName.c_str(), gmx::test::ExactTextMatch());

        gmx::test::StdioTestHelper stdioHelper(&fileManager());
        stdioHelper.redirectStringToStdin("q\n");

        ASSERT_EQ(0, gmx_make_ndx(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }
};

TEST_F(GmxMakeNdx, WritesDefaultProteinIndexGroups)
{
    std::string sysName("villin");
    runTest(sysName);
}

} // namespace
