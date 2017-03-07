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
 * Tests for gmx clustsize
 *
 * \todo These will be superseded by tests of the new style analysis
 * modules.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
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

class ClustsizeTest : public IntegrationTestFixture
{
    public:
        ClustsizeTest() : data_(), checker_(data_.rootChecker()) {}

        void runTest(const char *option, const char *filename)
        {
            CommandLine caller;
            caller.append("clustsize");

            caller.addOption("-f", fileManager_.getInputFilePath("clustsize.pdb"));
            caller.addOption("-n", fileManager_.getInputFilePath("clustsize.ndx"));

            std::string mcFileName = fileManager_.getTemporaryFilePath(filename);
            caller.addOption(option, mcFileName);

            EXPECT_EQ(0, gmx_clustsize(caller.argc(), caller.argv()));

            EXPECT_TRUE(Path::exists(mcFileName));
            TextInputFile    mcFile(mcFileName);
            XvgMatchSettings settings;
            auto             outputFilesChecker = checker_.checkCompound("OutputFiles", "Files");
            auto             fileChecker        = outputFilesChecker.checkCompound("File", "-o");
            checkXvgFile(&mcFile, &fileChecker, settings);
        }

        TestReferenceData    data_;
        TestReferenceChecker checker_;
};

TEST_F(ClustsizeTest, MaxClust)
{
    runTest("-mc", "maxclust.xvg");
}

TEST_F(ClustsizeTest, NClust)
{
    runTest("-nc", "nclust.xvg");
}

TEST_F(ClustsizeTest, AvClust)
{
    runTest("-ac", "avclust.xvg");
}

TEST_F(ClustsizeTest, HistoClust)
{
    runTest("-hc", "histoclust.xvg");
}


} // namespace

} // namespace

} // namespace
