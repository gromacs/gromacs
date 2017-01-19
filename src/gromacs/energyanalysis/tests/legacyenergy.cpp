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
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"

#include "testutils/cmdlinetest.h"
#include "testutils/integrationtests.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/xvgtest.h"

namespace gmx
{
namespace test
{
namespace
{

class DhdlTest : public IntegrationTestFixture
{
    public:
        DhdlTest() : data_(), checker_(data_.rootChecker()) {}

        void runTest()
        {
            CommandLine caller;
            caller.append("energy");

            caller.addOption("-s", fileManager_.getInputFilePath("dhdl.tpr"));
            caller.addOption("-f", fileManager_.getInputFilePath("dhdl.edr"));

            std::string xvgFileName = fileManager_.getTemporaryFilePath("dhdl.xvg");
            caller.addOption("-odh", xvgFileName);

            EXPECT_EQ(0, gmx_energy(caller.argc(), caller.argv()));

            EXPECT_TRUE(Path::exists(xvgFileName));
            TextInputFile    xvgFile(xvgFileName);
            XvgMatchSettings settings;
            auto             outputFilesChecker = checker_.checkCompound("OutputFiles", "Files");
            auto             fileChecker        = outputFilesChecker.checkCompound("File", "-o");
            checkXvgFile(&xvgFile, &fileChecker, settings);
        }

        TestReferenceData    data_;
        TestReferenceChecker checker_;
};

TEST_F(DhdlTest, ExtractDhdl)
{
    runTest();
}

class EnergyTest : public IntegrationTestFixture
{
    public:
        EnergyTest() : data_(), checker_(data_.rootChecker()) {}

        void runTest(const char *stringForStdin)
        {
            CommandLine caller;
            caller.append("energy");

            caller.addOption("-f", fileManager_.getInputFilePath("ener.edr"));

            std::string xvgFileName = fileManager_.getTemporaryFilePath("energy.xvg");
            caller.addOption("-o", xvgFileName);

            redirectStringToStdin(stringForStdin);
            EXPECT_EQ(0, gmx_energy(caller.argc(), caller.argv()));

            EXPECT_TRUE(Path::exists(xvgFileName));
            TextInputFile    xvgFile(xvgFileName);
            XvgMatchSettings settings;
            auto             outputFilesChecker = checker_.checkCompound("OutputFiles", "Files");
            auto             fileChecker        = outputFilesChecker.checkCompound("File", "-o");
            checkXvgFile(&xvgFile, &fileChecker, settings);
        }

        TestReferenceData    data_;
        TestReferenceChecker checker_;
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

class ViscosityTest : public IntegrationTestFixture
{
    public:
        ViscosityTest() : data_(), checker_(data_.rootChecker()) {}

        std::string addFileToTest(CommandLine *caller,
                                  const char  *optionName,
                                  const char  *outputFileName)
        {
            std::string xvgFileName = fileManager_.getTemporaryFilePath(outputFileName);
            caller->addOption(optionName, xvgFileName);
            return xvgFileName;
        }

        void runTest(const char *optionName,
                     const char *outputFileName)
        {
            CommandLine caller;
            caller.append("energy");

            caller.addOption("-f", fileManager_.getInputFilePath("ener.edr"));
            caller.addOption("-vis", fileManager_.getTemporaryFilePath("visco.xvg"));

            std::string xvgFileNameToTest;

            /* -vis can write a lot of non-conditional output files,
                so we use temporary paths to clean up files that are
                not the ones being tested in this test */
            if (0 == std::strcmp(optionName, "-evisco"))
            {
                xvgFileNameToTest = addFileToTest(&caller, optionName, outputFileName);
            }
            else
            {
                caller.addOption("-evisco", fileManager_.getTemporaryFilePath("evisco"));
            }
            if (0 == std::strcmp(optionName, "-eviscoi"))
            {
                xvgFileNameToTest = addFileToTest(&caller, optionName, outputFileName);
            }
            else
            {
                caller.addOption("-eviscoi", fileManager_.getTemporaryFilePath("eviscoi"));
            }
            if (0 == std::strcmp(optionName, "-corr"))
            {
                xvgFileNameToTest = addFileToTest(&caller, optionName, outputFileName);
            }
            else
            {
                caller.addOption("-corr", fileManager_.getTemporaryFilePath("enecorr"));
            }

            caller.addOption("-o", fileManager_.getTemporaryFilePath("energy"));

            EXPECT_EQ(0, gmx_energy(caller.argc(), caller.argv()));

            EXPECT_TRUE(Path::exists(xvgFileNameToTest));
            TextInputFile    xvgFile(xvgFileNameToTest);
            XvgMatchSettings settings;
            settings.tolerance = relativeToleranceAsFloatingPoint(1e-4, 1e-8);
            auto             outputFilesChecker = checker_.checkCompound("OutputFiles", "Files");
            auto             fileChecker        = outputFilesChecker.checkCompound("File", "-o");
            checkXvgFile(&xvgFile, &fileChecker, settings);
        }

        TestReferenceData    data_;
        TestReferenceChecker checker_;
};

TEST_F(ViscosityTest, EinsteinViscosity)
{
    runTest("-evisco", "evisco.xvg");
}

TEST_F(ViscosityTest, EinsteinViscosityIntegral)
{
    runTest("-eviscoi", "eviscoi.xvg");
}

} // namespace
} // namespace
} // namespace
