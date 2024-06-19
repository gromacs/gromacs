/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Tests for tabulated bonded interactions
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

//! Format string for building a configurable .top file
static const char* g_butaneTopFileFormatString =
        "\
[ defaults ]\n\
; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ\n\
  1		1		no		1.0	1.0\n\
\n\
[ atomtypes ]\n\
;name        mass      charge   ptype            c6           c12\n\
  CH2    14.02700       0.000       A   0.90975E-02   0.35333E-04\n\
  CH3    15.03500       0.000       A   0.88765E-02   0.26150E-04\n\
\n\
[ moleculetype ]\n\
;             name    nrexcl\n\
            butane   3\n\
\n\
[ atoms ]\n\
;   nr    type   resnr  residu    atom    cgnr\n\
     1     CH3       1     BUT      C1       1\n\
     2     CH2       1     BUT      C2       2\n\
     3     CH2       1     BUT      C3       3\n\
     4     CH3       1     BUT      C4       4\n\
\n\
%s\n\
\n\
[ system ]\n\
; The name of the system to be simulated\n\
A single butane\n\
\n\
[ molecules ]\n\
; Molname             Number\n\
Butane                   1\n\
";

//! Test fixture for bonded interactions
class BondedInteractionsTest : public gmx::test::MdrunTestFixture
{
public:
    //! Execute the trajectory writing test
    void setupGrompp(const char* interaction)
    {
        runner_.topFileName_ = fileManager_.getTemporaryFilePath("butane1.top").string();
        TextWriter::writeFileFromString(runner_.topFileName_,
                                        formatString(g_butaneTopFileFormatString, interaction));
        runner_.groFileName_ = gmx::test::TestFileManager::getInputFilePath("butane1.gro").string();
        runner_.ndxFileName_ = gmx::test::TestFileManager::getInputFilePath("butane1.ndx").string();
        runner_.useEmptyMdpFile();
    }
    //! Prepare an mdrun caller
    CommandLine setupMdrun()
    {
        CommandLine rerunCaller;
        rerunCaller.append("mdrun");
        rerunCaller.addOption("-rerun", runner_.groFileName_);
        return rerunCaller;
    }
    //! Check the output of mdrun
    void checkMdrun()
    {
        // TODO verifying some energies and forces would be good,
        // once other code in gerrit is reviewed
    }
};

// This test ensures that a normal non-tabulated bond interaction works
TEST_F(BondedInteractionsTest, NormalBondWorks)
{
    setupGrompp(
            "[ bonds ]\n\
;  ai    aj funct           c0           c1\n\
    1     2     1 1.530000e-01 3.347000e+05");
    EXPECT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupMdrun();
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    checkMdrun();
}

// This test ensures that a normal abulated bond interaction works
TEST_F(BondedInteractionsTest, TabulatedBondWorks)
{
    setupGrompp(
            "[ bonds ]\n\
;  ai    aj funct  n     k\n\
    1     2     8  0  1000");
    EXPECT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupMdrun();
    std::string tableFileName = gmx::test::TestFileManager::getInputFilePath("butane_b0.xvg").string();
    rerunCaller.addOption("-tableb", tableFileName);
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    checkMdrun();
}

// This test ensures that a normal non-tabulated angle interaction works
TEST_F(BondedInteractionsTest, NormalAngleWorks)
{
    setupGrompp(
            "[ angles ]\n\
;  ai    aj    ak funct           c0           c1\n\
    1     2     3     1 1.110000e+02 4.602000e+02");
    EXPECT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupMdrun();
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    checkMdrun();
}

// This test ensures that a tabulated angle interaction works
TEST_F(BondedInteractionsTest, TabulatedAngleWorks)
{
    setupGrompp(
            "[ angles ]\n\
;  ai    aj    ak funct  n     k\n\
    1     2     3     8  0  1000");
    EXPECT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupMdrun();
    std::string tableFileName = gmx::test::TestFileManager::getInputFilePath("butane_a0.xvg").string();
    rerunCaller.addOption("-tableb", tableFileName);
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    checkMdrun();
}

// This test ensures that a normal non-tabulated dihedral interaction works
TEST_F(BondedInteractionsTest, NormalDihedralWorks)
{
    setupGrompp(
            "[ dihedrals ]\n \
;  ai    aj    ak    al funct     c0     c1     c2      c3     c4      c5\n\
    1     2     3     4     3 9.2789 12.156 -13.12 -3.0597  26.24 -31.495");
    EXPECT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupMdrun();
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    checkMdrun();
}

// This test ensures that a tabulated dihedral interaction works
TEST_F(BondedInteractionsTest, TabulatedDihedralWorks)
{
    setupGrompp(
            "[ dihedrals ]\n\
;  ai    aj    ak    al funct   n     k\n\
    1     2     3     4     8   0  1000");
    EXPECT_EQ(0, runner_.callGrompp());

    test::CommandLine rerunCaller = setupMdrun();
    std::string tableFileName = gmx::test::TestFileManager::getInputFilePath("butane_d0.xvg").string();
    rerunCaller.addOption("-tableb", tableFileName);
    ASSERT_EQ(0, runner_.callMdrun(rerunCaller));
    checkMdrun();
}

} // namespace test

} // namespace gmx
