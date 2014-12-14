/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Tests for analysis of energy files
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "gromacs/energyanalysis/energyinfo.h"

//#include "moduletest.h"

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"

#include "../dhdl.h"
#include "../fluctprops.h"
#include "../freeenergydifference.h"
#include "../simple.h"
#include "../viscosity.h"

namespace
{

using gmx::test::CommandLine;

class SimpleTest : public gmx::test::CommandLineTestBase
{
    public:
        SimpleTest()
        {
            setInputFile("-f", "ener.edr");
            setOutputFile("-o", "energy.xvg");
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            gmx::test::TestReferenceChecker rootChecker(this->rootChecker());
            rootChecker.checkString(args.toString(), "CommandLine");
            int rc = gmx::test::CommandLineTestHelper::runModule(
                        &gmx::SimpleInfo::create, &cmdline);
            ASSERT_EQ(0, rc);
            checkOutputFiles();
        }
};

TEST_F(SimpleTest, ExtractsEnergy)
{
    const char *const cmdline[] = {
        "energy", "-xvg", "none", "-term", "'Potential Kinetic-En. Total-Energy'"
    };
    //setInputFile("-f", "ener.edr");
    //  setOutputFile("-o", "energy.xvg");
    runTest(CommandLine(cmdline));
}

//! Test fixture for the simple energyanalysis module.
/*typedef gmx::test::EnergyAnalysisModuleTestFixture<gmx::SimpleInfo>
    SimpleTest2;

   TEST_F(SimpleTest2, ExtractsEnergy2)
   {
    const char *const cmdline[] = {
        "energy", "-term", "Potential Kinetic-En. Total-Energy"
    };
    setEnergy("ener.edr");
    setOutputFile("-o", "energy.xvg");
    runTest(CommandLine(cmdline));
   }
 */


class FluctPropsTest : public gmx::test::CommandLineTestBase
{
    public:
        FluctPropsTest()
        {
            setInputFile("-f", "ener.edr");
            setOutputFile("-convergence", "convergence.xvg");
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            gmx::test::TestReferenceChecker rootChecker(this->rootChecker());
            rootChecker.checkString(args.toString(), "CommandLine");
            ASSERT_EQ(0, gmx::test::CommandLineTestHelper::runModule(
                              &gmx::FluctPropsInfo::create, &cmdline));

            checkOutputFiles();
        }
};

TEST_F(FluctPropsTest, ExtractFluctuations)
{
    const char *const cmdline[] = {
        "fluctprops", "-xvg", "none", "-nmol", "1000"
    };
    runTest(CommandLine(cmdline));
}

class ViscosityTest : public gmx::test::CommandLineTestBase
{
    public:
        ViscosityTest()
        {
            setInputFile("-f", "ener.edr");
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            gmx::test::TestReferenceChecker rootChecker(this->rootChecker());
            rootChecker.checkString(args.toString(), "CommandLine");
            ASSERT_EQ(0, gmx::test::CommandLineTestHelper::runModule(
                              &gmx::ViscosityInfo::create, &cmdline));

            checkOutputFiles();
        }
};

TEST_F(ViscosityTest, Einstein)
{
    const char *const cmdline[] = {
        "viscosity", "-xvg", "none"
    };
    setOutputFile("-evis", "einstein.xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(ViscosityTest, EinsteinIntegral)
{
    const char *const cmdline[] = {
        "viscosity", "-xvg", "none"
    };
    setOutputFile("-eivis", "einsteinintegral.xvg");
    runTest(CommandLine(cmdline));
}

TEST_F(ViscosityTest, Fanourgakis)
{
    const char *const cmdline[] = {
        "viscosity", "-xvg", "none", "-pcorr", "presscorr"
    };
    runTest(CommandLine(cmdline));
}

} // namespace
