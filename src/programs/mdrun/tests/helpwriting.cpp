/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * This implements tests on mdrun help writing.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"

#include "programs/mdrun/mdrun_main.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(MdrunTest, WritesHelp)
{
    // Make a stream to which we want gmx mdrun -h to write the help.
    StringOutputStream outputStream;
    TextWriter         writer(&outputStream);

    // Use that stream to set up a global help context. Legacy tools
    // like mdrun call parse_common_args, which recognizes the
    // existence of a global help context. That context triggers the
    // writing of help and a fast exit of the tool.
    HelpLinks*                   links = nullptr;
    CommandLineHelpContext       context(&writer, eHelpOutputFormat_Console, links, "dummy");
    GlobalCommandLineHelpContext global(context);

    // Call mdrun to get the help printed to the stream
    CommandLine caller;
    caller.append("mdrun");
    caller.append("-h");
    gmx_mdrun(caller.argc(), caller.argv());

    // Check whether the stream matches the reference copy.
    TestReferenceData    refData;
    TestReferenceChecker checker(refData.rootChecker());
    checker.checkString(outputStream.toString(), "Help string");
};

} // namespace
} // namespace test
} // namespace gmx
