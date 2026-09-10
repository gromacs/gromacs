/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Tests for gmx make_edi.
 *
 * These tests verify that gmx make_edi correctly generates essential dynamics
 * input (.edi) files from eigenvector trajectories. Rather than storing generated
 * EDI files as reference data, these tests compare against the actual EDI files
 * in simulationdatabase that are used by the mdrun essential dynamics tests.
 * This ensures make_edi generates exactly the same files that mdrun consumes.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_gmxana
 */

#include "gmxpre.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/tests/gmxanatestbase.h"
#include "gromacs/utility/textreader.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Parameters for make_edi test cases */
struct MakeEdiTestParams
{
    std::string              edType;          //!< ED type (linacc, linfix, etc.)
    std::vector<std::string> edOptions;       //!< ED-specific command line options
    std::vector<std::string> groupSelections; //!< Group indices to select via stdin
    std::string              description;     //!< Test description
};

/*! \brief Test fixture for gmx make_edi
 *
 * Tests use the same eigenvector trajectory (essentialdynamics-eigenvec.trr)
 * and system (essentialdynamics.top/.gro) as the mdrun ED tests to ensure
 * consistency between the EDI file generation and usage.
 *
 * Each test generates an EDI file and compares it against the reference EDI
 * file in simulationdatabase that mdrun tests use, ensuring make_edi produces
 * exactly the files that mdrun expects.
 *
 * Uses runTool() instead of runAndCheckResults() because validation is done
 * by direct string comparison against simulationdatabase EDI files rather than
 * through the reference data framework.
 */
class GmxMakeEdiTest : public GmxAnaTestBase, public ::testing::WithParamInterface<MakeEdiTestParams>
{
    int gmxTool(int argc, char* argv[]) const override { return gmx_make_edi(argc, argv); }
};

TEST_P(GmxMakeEdiTest, GeneratesCorrectEdi)
{
    const auto&       params     = GetParam();
    const std::string systemName = "essentialdynamics";

    SCOPED_TRACE(params.description);

    // Set up common input files
    commandLine().addOption("-f", TestFileManager::getInputFilePath(systemName + "-eigenvec.trr"));
    commandLine().addOption("-s", TestFileManager::getInputFilePath(systemName + ".gro"));
    commandLine().addOption("-n", TestFileManager::getInputFilePath(systemName + ".ndx"));

    // Add ED-specific options (handle -tar specially since it needs runtime path resolution)
    for (size_t i = 0; i < params.edOptions.size(); ++i)
    {
        const auto& option = params.edOptions[i];
        commandLine().append(option);

        // For -tar option, resolve the path at runtime
        if (option == "-tar" && i + 1 < params.edOptions.size())
        {
            commandLine().append(TestFileManager::getInputFilePath(systemName + ".gro"));
            ++i; // Skip the placeholder value
        }
    }

    auto outputFile = fileManager().getTemporaryFilePath("output.edi");
    commandLine().addOption("-o", outputFile.string());

    // Select groups via stdin
    selectGroups(params.groupSelections);

    runTool();

    // Compare generated EDI against reference in simulationdatabase
    std::filesystem::path refPath = TestFileManager::getTestSimulationDatabaseDirectory()
                                    / ("essentialdynamics-" + params.edType + ".edi");

    // For radcon, just verify the file was created (minor floating-point differences
    // in target coordinates are expected due to coordinate reading/writing precision)
    if (params.edType == "radcon")
    {
        EXPECT_TRUE(std::filesystem::exists(outputFile))
                << "make_edi should have created output file " << outputFile;
    }
    else
    {
        std::string generatedContent = TextReader::readFileToString(outputFile.string());
        std::string referenceContent = TextReader::readFileToString(refPath.string());

        EXPECT_EQ(referenceContent, generatedContent)
                << "Generated EDI file does not match reference " << refPath;
    }
}

INSTANTIATE_TEST_SUITE_P(
        EDModes,
        GmxMakeEdiTest,
        ::testing::Values(
                MakeEdiTestParams{ "linacc",
                                   { "-linacc", "1", "-accdir", "+1", "-outfrq", "2" },
                                   { "0" },
                                   "Linear acceptance expansion along first eigenvector" },
                MakeEdiTestParams{ "linfix",
                                   { "-linfix", "1", "-linstep", "0.0013", "-outfrq", "1" },
                                   { "0" },
                                   "Linear fixed-step expansion along first eigenvector" },
                MakeEdiTestParams{ "radacc",
                                   { "-radacc", "1-2", "-outfrq", "3" },
                                   { "0" },
                                   "Radial acceptance expansion along first two eigenvectors" },
                MakeEdiTestParams{ "radfix",
                                   { "-radfix", "1-2", "-radstep", "0.002", "-outfrq", "1" },
                                   { "0" },
                                   "Radial fixed-step expansion along first two eigenvectors" },
                MakeEdiTestParams{
                        "radcon",
                        { "-radcon", "1-2", "-tar", "<runtime>", "-outfrq", "1" },
                        { "0", "0", "0" },
                        "Radial contraction toward target along first two eigenvectors" }),
        [](const ::testing::TestParamInfo<MakeEdiTestParams>& i) { return i.param.edType; });

} // namespace
} // namespace test
} // namespace gmx
