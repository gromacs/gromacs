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
 * Tests for solvation.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "../solvate.h"
#include "testutils/integrationtests.h"
#include "testutils/cmdlinetest.h"
#include "gromacs/fileio/futil.h"

namespace
{

//! Helper typedef
typedef int (*CMainFunction)(int argc, char *argv[]);

//! Test fixture for gmx solvate
class SolvateTest : public gmx::test::IntegrationTestFixture
{
    public:
        //! Constructor
        SolvateTest()
            : cpFileName(fileManager_.getInputFilePath("spc-and-methanol.gro")),
              topFileName(fileManager_.getInputFilePath("spc-and-methanol.top")),
              outputFileName(fileManager_.getTemporaryFilePath("out.gro")),
              theFunction(gmx_solvate)
        {
            caller.append("solvate");
            caller.addOption("-o", outputFileName);
        }

    public:
        //! Name of file to use for -cp
        std::string            cpFileName;
        //! Name of input file to use for -p (if used)
        std::string            topFileName;
        //! Name of output file to use for -o
        std::string            outputFileName;
        //! Helper object for managing the call to gmx_solvate
        gmx::test::CommandLine caller;
        //! Points to the function to be tested
        CMainFunction          theFunction;
};

TEST_F(SolvateTest, cs_box_Works)
{
    caller.append("-cs"); // use default solvent box
    caller.addOption("-box", "1.1");

    ASSERT_EQ(0, theFunction(caller.argc(), caller.argv()));
}

TEST_F(SolvateTest, cs_cp_Works)
{
    caller.append("-cs"); // use default solvent box
    caller.addOption("-cp", cpFileName);

    ASSERT_EQ(0, theFunction(caller.argc(), caller.argv()));
}

TEST_F(SolvateTest, cs_cp_p_Works)
{
    caller.append("-cs"); // use default solvent box
    caller.addOption("-cp", cpFileName);

    std::string modifiableTopFileName = fileManager_.getTemporaryFilePath(".top");
    gmx_file_copy(topFileName.c_str(), modifiableTopFileName.c_str(), true);
    caller.addOption("-p", modifiableTopFileName);

    ASSERT_EQ(0, theFunction(caller.argc(), caller.argv()));
}

} // namespace
