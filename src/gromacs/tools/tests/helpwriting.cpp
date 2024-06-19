/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * This implements tests on tool help writing. Based on mdrun test version.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/tools/convert_tpr.h"
#include "gromacs/tools/dump.h"
#include "gromacs/tools/report_methods.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"

namespace gmx
{
namespace test
{
namespace
{

class HelpwritingTest : public gmx::test::CommandLineTestBase
{
public:
    void runTest(gmx::ICommandLineModule* module) { testWriteHelp(module); }
};

TEST_F(HelpwritingTest, ConvertTprWritesHelp)
{
    const std::unique_ptr<gmx::ICommandLineModule> module(gmx::ICommandLineOptionsModule::createModule(
            "convert-tpr", "Dummy Info", ConvertTprInfo::create()));
    runTest(module.get());
};


TEST_F(HelpwritingTest, DumpWritesHelp)
{
    const std::unique_ptr<gmx::ICommandLineModule> module(
            gmx::ICommandLineOptionsModule::createModule("dump", "Dummy Info", DumpInfo::create()));
    runTest(module.get());
};

TEST_F(HelpwritingTest, ReportMethodsWritesHelp)
{
    const std::unique_ptr<gmx::ICommandLineModule> module(gmx::ICommandLineOptionsModule::createModule(
            "report-methods", "Dummy Info", ReportMethodsInfo::create()));
    runTest(module.get());
};

} // namespace
} // namespace test
} // namespace gmx
