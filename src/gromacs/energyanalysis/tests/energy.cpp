/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "gromacs/energyanalysis/modules/energy.h"

#include <cstdlib>

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/energyanalysis/energyanalysisrunner.h"
#include "gromacs/energyanalysis/modules/dhdl.h"
#include "gromacs/energyanalysis/modules/fluctprops.h"
#include "gromacs/energyanalysis/modules/viscosity.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/snprintf.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

namespace gmx
{

namespace energyanalysis
{

namespace
{

using test::CommandLine;

class EnergyTest : public gmx::test::CommandLineTestBase
{
    public:
        EnergyTest()
        {
            setInputFile("-f", "ener.edr");
        }

        /*! \brief Run the EnergyTest
         *
         * \param[in] args      The command line arguments
         * \param[in] option    Command line flag
         * \param[in] filename  The output file name
         */
        void runTest(const CommandLine &args,
                     const char        *option,
                     const char        *filename)
        {
            double            tolerance = 1e-4;
            test::XvgMatch    xvg;
            setOutputFile(option, filename,
                          xvg.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance)));
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            rootChecker().checkString(args.toString(), "CommandLine");
            ICommandLineOptionsModulePointer runner(EnergyAnalysisRunner::createModule(EnergyInfo::create()));

            ASSERT_EQ(0, gmx::test::CommandLineTestHelper::runModuleDirect(std::move(runner), &cmdline));

            checkOutputFiles();
        }
};

TEST_F(EnergyTest, ExtractEnergy)
{
    const char *const cmdline[] = {
        "energy", "-term", "Potential Kinetic-En. Total-Energy"
    };
    runTest(CommandLine(cmdline), "-o", "energy.xvg");
}

TEST_F(EnergyTest, ExtractEnergyByNumber)
{
    const char *const cmdline[] = {
        "energy", "-term", "4 6 9"
    };
    runTest(CommandLine(cmdline), "-o", "energy.xvg");
}

TEST_F(EnergyTest, ExtractEnergyMixed)
{
    const char *const cmdline[] = {
        "energy", "-term", "Pressu 7 box-z vol"
    };
    runTest(CommandLine(cmdline), "-o", "energy.xvg");
}

} // namespace

} // namespace energyanalysis

} // namespace gmx
