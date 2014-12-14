/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/energyanalysis/modules/viscosity.h"

#include <cstdlib>

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/energyanalysis/energyanalysisrunner.h"

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

class ViscosityTest : public gmx::test::CommandLineTestBase
{
    public:
        ViscosityTest()
        {
            setInputFile("-f", "ener.edr");
        }

        /*! \brief Run the ViscosityTest
         *
         * \param[in] option      The command line option
         * \param[in] filename    The corresponding output file
         */
        void runTest(const char *option, const char *filename, bool noUsePDiff)
        {
            const char *const command[] = { "viscosity" };
            double            tolerance = 1e-4;
            test::XvgMatch    xvg;
            setOutputFile(option, filename,
                          xvg.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance)));
            CommandLine &cmdline = commandLine();
            CommandLine  args    = CommandLine(command);
            cmdline.merge(args);
            if (noUsePDiff)
            {
                const char *const extra_args[]   = { "-nousePressureTraceDiff" };
                CommandLine       extra          = CommandLine(extra_args);
                cmdline.merge(extra);
            }
            rootChecker().checkString(args.toString(), "CommandLine");
            ICommandLineOptionsModulePointer runner(EnergyAnalysisRunner::createModule(ViscosityInfo::create()));

            ASSERT_EQ(0, gmx::test::CommandLineTestHelper::runModuleDirect(std::move(runner), &cmdline));

            checkOutputFiles();
        }
};

TEST_F(ViscosityTest, EinsteinViscosity)
{
    runTest("-evisco", "einstein.xvg", true);
}

TEST_F(ViscosityTest, EinsteinViscosityAllPressureElements)
{
    runTest("-evisco", "einstein.xvg", false);
}

TEST_F(ViscosityTest, EinsteinViscosityIntegral)
{
    runTest("-eviscoi", "einsteinintegral.xvg", true);
}

TEST_F(ViscosityTest, EinsteinViscosityIntegralAllPressureElements)
{
    runTest("-eviscoi", "einsteinintegral.xvg", false);
}

}

} // namespace energyanalysis

} // namespace gmx
