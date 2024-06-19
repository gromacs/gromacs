/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests for functionality of the "extract-cluster" trajectory analysis module.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/extract_cluster.h"

#include <array>
#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

#include "moduletest.h"

namespace gmx
{

namespace test
{

namespace
{

/********************************************************************
 * Tests for gmx::analysismodules::ExtractCluster.
 */

//! Helper struct to combine filename, path and FileMatcher
struct ManualOutputFile
{
    //! Generated file name.
    std::string filename;
    //! Full generated file path.
    std::string fullFilepath;
    //! Corresponding file matcher to compare to reference.
    FileMatcherPointer matcher;
};

//! Test fixture for the convert-trj analysis module.
class ExtractClusterModuleTest : public AbstractTrajectoryAnalysisModuleTestFixture
{
protected:
    TrajectoryAnalysisModulePointer createModule() override
    {
        return analysismodules::ExtractClusterInfo::create();
    }

public:
    //! Constructor
    ExtractClusterModuleTest();

    //! Compare generated files to expected ones.
    void compareFiles();

private:
    //! Array of names and paths for generated files.
    std::array<ManualOutputFile, 8> generatedFiles;
};

ExtractClusterModuleTest::ExtractClusterModuleTest()
{
    auto outputFilepath = std::filesystem::current_path();

    fileManager().setOutputTempDirectory(outputFilepath);
    // Those are for cleaning up the files generated during testing.
    fileManager().getTemporaryFilePath("csizew.xpm");
    fileManager().getTemporaryFilePath("csize.xpm");

    setTrajectory("extract_cluster.trr");
    setInputFile("-clusters", "extract_cluster.ndx");

    int fileNumber = 1;
    for (auto& generatedFile : generatedFiles)
    {
        generatedFile.filename = gmx::concatenateBeforeExtension(
                                         "test.g96", gmx::formatString("_Cluster_000%d", fileNumber))
                                         .string();
        generatedFile.matcher      = TextFileMatch(ExactTextMatch()).createFileMatcher();
        generatedFile.fullFilepath = fileManager().getTemporaryFilePath(generatedFile.filename).string();
        fileNumber++;
    }
}

void ExtractClusterModuleTest::compareFiles()
{
    TestReferenceChecker outputChecker(rootChecker().checkCompound("ClusterOutputFiles", "Files"));
    for (const auto& file : generatedFiles)
    {
        TestReferenceChecker manualFileChecker(outputChecker.checkCompound("File", file.filename.c_str()));
        file.matcher->checkFile(file.fullFilepath, &manualFileChecker);
    }
}

TEST_F(ExtractClusterModuleTest, WorksWithAllAtoms)
{
    std::string       realFileName = TestFileManager::getTestSpecificFileName("test.g96").string();
    const char* const cmdline[]    = { "extract-cluster", "-o", realFileName.c_str() };

    runTest(CommandLine(cmdline));
    compareFiles();
}

TEST_F(ExtractClusterModuleTest, WorksWithAtomSubset)
{
    std::string       realFileName = TestFileManager::getTestSpecificFileName("test.g96").string();
    const char* const cmdline[]    = {
        "extract-cluster", "-o", realFileName.c_str(), "-select", "atomnr 1 2"
    };

    runTest(CommandLine(cmdline));
    compareFiles();
}

} // namespace
} // namespace test
} // namespace gmx
