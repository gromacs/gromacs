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

#include <cstdlib>

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/energyanalysis/dhdl.h"
#include "gromacs/energyanalysis/energyinfo.h"
#include "gromacs/energyanalysis/fluctprops.h"
#include "gromacs/energyanalysis/freeenergydifference.h"
#include "gromacs/energyanalysis/simple.h"
#include "gromacs/energyanalysis/viscosity.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/snprintf.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/xvgtests.h"

namespace gmx
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
         * \param[in] tolerance The relative tolerance in the output
         * \param[in] module    The module to test
         */
        void runTest(const CommandLine        &args,
                     const double              tolerance,
                     ICommandLineOptionsModule *(*module)())
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            gmx::test::TestReferenceChecker rootChecker(this->rootChecker());
            rootChecker.checkString(args.toString(), "CommandLine");
            rootChecker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
            ASSERT_EQ(0, gmx::test::CommandLineTestHelper::runModule(
                              module, &cmdline));
            checkOutputFiles();
        }
};

TEST_F(EnergyTest, SimpleExtractEnergy)
{
    const char *const cmdline[] = {
        "energy", "-term", "'Potential Kinetic-En. Total-Energy'"
    };
    double            tolerance = 0.001;
    setOutputFileXvg("-o", "energy.xvg");
    runTest(CommandLine(cmdline), tolerance, &SimpleInfo::create);
}

TEST_F(EnergyTest, ExtractFluctuations)
{
    const char *const cmdline[] = { "fluctprops", "-nmol", "1000" };
    double            tolerance = 0.001;
    setOutputFileXvg("-convergence", "convergence.xvg");
    runTest(CommandLine(cmdline), tolerance, &FluctPropsInfo::create);
}

TEST_F(EnergyTest, EinsteinViscosity)
{
    const char *const cmdline[] = { "viscosity" };
    double            tolerance = 0.001;
    setOutputFileXvg("-evis", "einstein.xvg");
    runTest(CommandLine(cmdline), tolerance, &ViscosityInfo::create);
}

TEST_F(EnergyTest, EinsteinViscosityIntegral)
{
    const char *const cmdline[] = { "viscosity" };
    double            tolerance = 0.001;
    setOutputFileXvg("-eivis", "einsteinintegral.xvg");
    runTest(CommandLine(cmdline), tolerance, &ViscosityInfo::create);
}

TEST_F(EnergyTest, FanourgakisViscosity)
{
    const char *const cmdline[] = { "viscosity" };
    double            tolerance = 0.001;
    setOutputFileXvg("-pcorr", "presscorr.xvg");
    runTest(CommandLine(cmdline), tolerance, &ViscosityInfo::create);
}

} // namespace

} // namespace gmx
