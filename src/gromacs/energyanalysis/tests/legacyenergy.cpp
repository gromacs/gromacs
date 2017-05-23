/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for gmx energy
 *
 * \todo These will be superseded by tests of the energyanalysis
 * modules.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include <cstring>

#include <string>

#include "gromacs/gmxana/gmx_ana.h"

#include "testutils/cmdlinetest.h"
#include "testutils/stdiohelper.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

namespace gmx
{
namespace test
{
namespace
{

class DhdlTest : public CommandLineTestBase
{
    public:
        void runTest()
        {
            auto &cmdline = commandLine();
            cmdline.append("energy");

            setInputFile("-s", "dhdl.tpr");
            setInputFile("-f", "dhdl.edr");
            setOutputFile("-odh", "dhdl.xvg", XvgMatch());

            ASSERT_EQ(0, gmx_energy(cmdline.argc(), cmdline.argv()));

            checkOutputFiles();
        }
};

TEST_F(DhdlTest, ExtractDhdl)
{
    runTest();
}

class OriresTest : public CommandLineTestBase
{
    public:
        void runTest(const char *stringForStdin)
        {
            auto &cmdline = commandLine();
            cmdline.append("energy");

            setInputFile("-s", "orires.tpr");
            setInputFile("-f", "orires.edr");
            test::XvgMatch    xvg;
            test::XvgMatch   &toler     = xvg.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, 1e-4));

            setOutputFile("-oten", ".xvg", toler);
            setOutputFile("-ora", ".xvg", toler);
            setOutputFile("-ort", ".xvg", toler);
            setOutputFile("-oda", ".xvg", toler);
            setOutputFile("-odr", ".xvg", toler);

            StdioTestHelper stdioHelper(&fileManager());
            stdioHelper.redirectStringToStdin(stringForStdin);
            ASSERT_EQ(0, gmx_nmr(cmdline.argc(), cmdline.argv()));

            checkOutputFiles();
        }
};

TEST_F(OriresTest, ExtractOrires)
{
    runTest("-1\n");
}

class EnergyTest : public CommandLineTestBase
{
    public:
        void runTest(const char *stringForStdin)
        {
            auto &cmdline = commandLine();
            cmdline.append("energy");

            setInputFile("-f", "ener.edr");
            setOutputFile("-o", "energy.xvg", XvgMatch());

            StdioTestHelper stdioHelper(&fileManager());
            stdioHelper.redirectStringToStdin(stringForStdin);
            ASSERT_EQ(0, gmx_energy(cmdline.argc(), cmdline.argv()));

            checkOutputFiles();
        }
};

TEST_F(EnergyTest, ExtractEnergy)
{
    runTest("Potential\nKinetic-En.\nTotal-Energy\n");
}

TEST_F(EnergyTest, ExtractEnergyByNumber)
{
    runTest("4 6 9");
}

TEST_F(EnergyTest, ExtractEnergyMixed)
{
    runTest("Pressu\n7\nbox-z\nvol\n");
}

class ViscosityTest : public CommandLineTestBase
{
    public:
        void runTest()
        {
            auto &cmdline = commandLine();
            setInputFile("-f", "ener.edr");
            setOutputFile("-vis", "visco.xvg", NoTextMatch());
            setOutputFile("-o", "energy.xvg", NoTextMatch());

            /* -vis can write a lot of non-conditional output files,
                so we use temporary paths to clean up files that are
                not the ones being tested in this test */
            if (!cmdline.contains("-evisco"))
            {
                setOutputFile("-evisco", "evisco.xvg", NoTextMatch());
            }
            if (!cmdline.contains("-eviscoi"))
            {
                setOutputFile("-eviscoi", "eviscoi.xvg", NoTextMatch());
            }
            if (!cmdline.contains("-corr"))
            {
                setOutputFile("-corr", "corr.xvg", NoTextMatch());
            }

            ASSERT_EQ(0, gmx_energy(cmdline.argc(), cmdline.argv()));

            checkOutputFiles();
        }
};

TEST_F(ViscosityTest, EinsteinViscosity)
{
    auto tolerance = relativeToleranceAsFloatingPoint(1e-4, 1e-5);
    setOutputFile("-evisco", "evisco.xvg", XvgMatch().tolerance(tolerance));
    runTest();
}

TEST_F(ViscosityTest, EinsteinViscosityIntegral)
{
    auto tolerance = relativeToleranceAsFloatingPoint(1e-4, 1e-5);
    setOutputFile("-eviscoi", "eviscoi.xvg", XvgMatch().tolerance(tolerance));
    runTest();
}

} // namespace
} // namespace test
} // namespace gmx
