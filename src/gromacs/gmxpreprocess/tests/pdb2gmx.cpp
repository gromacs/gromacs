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
 * Tests for pdb2gmx
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/pdb2gmx.h"

#include "gromacs/utility/futil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::ConfMatch;

class Pdb2gmxTest : public gmx::test::CommandLineTestBase
{
    public:
        Pdb2gmxTest()
        {
            // TODO It would be preferable to just scrub the content
            // that actually varies, but we don't have enough regular
            // expression support for that yet.
            textMatcher_.addLineToSkip("^;[[:blank:]] *File '.*' was generated.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*By user:.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*On host:.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*At date:.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*:-\\).*\\(-:.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*Executable:.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*Data prefix:.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*Working dir:.*\n");
            textMatcher_.addLineToSkip("^;[[:blank:]]*gmxpreprocess-test.*\n");
            setOutputFile("-o", "conf.gro", ConfMatch());
            setOutputFile("-p", "topol.top", TextFileMatch(textMatcher_));
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            ASSERT_EQ(0, gmx_pdb2gmx(cmdline.argc(), cmdline.argv()));
            checkOutputFiles();
        }
        FilteringExactTextMatch textMatcher_;
};

TEST_F(Pdb2gmxTest, BasicWorks)
{
    const char *const cmdline[] = {
        "pdb2gmx", "-ignh", "-ff", "oplsaa", "-water", "tip3p", "-vsite", "none"
    };
    setInputFile("-f", "1aml.pdb");
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace
} // namespace
