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
 * Tests for insertion of molecules.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "../insert-molecules.h"
#include "testutils/integrationtests.h"
#include "testutils/cmdlinetest.h"
#include "gromacs/fileio/futil.h"

namespace
{

//! Helper typedef
typedef int (*CMainFunction)(int argc, char *argv[]);

//! Test fixture for insert-molecules
class InsertMoleculesTest : public gmx::test::IntegrationTestFixture
{
    public:
        //! Constructor
        InsertMoleculesTest()
            : ciFileName(fileManager_.getInputFilePath("x.gro")),
              outputFileName(fileManager_.getTemporaryFilePath("out.gro")),
              theFunction(gmx_insert_molecules)
        {
            caller.append("insert-molecules");
            caller.addOption("-o", outputFileName);
        }

    public:
        //! Name of file to use for -ci
        std::string            ciFileName;
        //! Name of output file to use for -o
        std::string            outputFileName;
        //! Helper object for managing the call to gmx_insert_molecules
        gmx::test::CommandLine caller;
        //! Points to the function to be tested
        CMainFunction          theFunction;
};

TEST_F(InsertMoleculesTest, f_ci_Works)
{
    caller.addOption("-f", fileManager_.getInputFilePath("spc-and-methanol.gro"));
    caller.addOption("-nmol", "1");
    caller.addOption("-ci", ciFileName);

    ASSERT_EQ(0, theFunction(caller.argc(), caller.argv()));
}

TEST_F(InsertMoleculesTest, box_ci_Works)
{
    caller.addOption("-box", "4");
    caller.addOption("-nmol", "5");
    caller.addOption("-ci", ciFileName);

    ASSERT_EQ(0, theFunction(caller.argc(), caller.argv()));
}

TEST_F(InsertMoleculesTest, f_box_ci_Works)
{
    caller.addOption("-f", fileManager_.getInputFilePath("spc-and-methanol.gro"));
    caller.addOption("-box", "4");
    caller.addOption("-nmol", "2");
    caller.addOption("-ci", ciFileName);

    ASSERT_EQ(0, theFunction(caller.argc(), caller.argv()));
}

// TODO Someone who knows what -ip is good for should write something
// to test it

} // namespace
