/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for functionality of the NNPotTopologyPreprocessor.
 * Adapted from qmmm/tests/qmmmtopologypreprocessor.cpp
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "gromacs/applied_forces/nnpot/nnpottopologypreprocessor.h"

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/cmdlinetest.h"
#include "testutils/loggertest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

enum class PbcType : int;

namespace gmx
{

namespace test
{

/*! \brief Test fixture for NNPotTopologyPreprocessor tests
 *
 * Provides a temporary file manager and a gmx_mtop_t struct for testing.
 */
class NNPotTopologyPreprocessorTest : public ::testing::Test
{
public:
    /*! \brief Generates tpr file from *.top and *.gro existing in the simulation database directory
     * and custom \c mdpContent and loads gmx_mtop_t from it.
     *
     * \param[in] simulationName Name of the simulation database directory
     * \param[in] mdpContent     Content of the mdp file
     */
    static std::unique_ptr<gmx_mtop_t> makeMtopFromFile(const std::string& simulationName,
                                                        const std::string& mdpContent)
    {
        const std::filesystem::path simData =
                gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();
        TestFileManager fileManager;

        // Generate empty mdp file
        const std::string mdpInputFileName =
                fileManager.getTemporaryFilePath(simulationName + ".mdp").string();
        gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpContent);

        // Generate tpr file
        const std::string tprName = fileManager.getTemporaryFilePath(simulationName + ".tpr").string();
        {
            gmx::test::CommandLine caller;
            caller.append("grompp");
            caller.addOption("-f", mdpInputFileName);
            caller.addOption("-p", (simData / simulationName).replace_extension(".top").string());
            caller.addOption("-c", (simData / simulationName).replace_extension(".gro").string());
            caller.addOption("-o", tprName);
            EXPECT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
        }

        // Load topology
        bool                        fullTopology;
        PbcType                     pbcType;
        matrix                      box;
        std::unique_ptr<gmx_mtop_t> mtop(std::make_unique<gmx_mtop_t>());
        readConfAndTopology(tprName.c_str(), &fullTopology, mtop.get(), &pbcType, nullptr, nullptr, box);
        return mtop;
    }

protected:
    LoggerTestHelper logHelper_;
};

TEST_F(NNPotTopologyPreprocessorTest, CanConstruct)
{
    std::vector<Index> nnpIndices = { 0, 1, 2 };
    EXPECT_NO_THROW(NNPotTopologyPreprocessor topPrep(nnpIndices));
}

TEST_F(NNPotTopologyPreprocessorTest, FourWatersFirstInQMRegion)
{
    // Reference input 4x SPCE waters from database 4waters.top
    // First water is NNP input
    std::vector<Index>          nnpAtomIndices = { 0, 1, 2 };
    std::unique_ptr<gmx_mtop_t> mtop           = makeMtopFromFile("4water", "");

    NNPotTopologyPreprocessor topPrep(nnpAtomIndices);
    MDLogger                  logger;
    WarningHandler            wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Neural network potential Interface is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded NNP atoms: 3\nNumber of regular atoms: 9\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 3\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of settles removed: 1 \\(replaced by 2 "
                                        "InteractionFunction::ConnectBonds\\) \n");
    topPrep.preprocess(mtop.get(), logHelper_.logger(), &wi);
}

TEST_F(NNPotTopologyPreprocessorTest, FourWatersSecondAndFourthInQMRegion)
{
    // Reference input 4x SPCE waters from database 4waters.top
    // second and fourth are NNP input
    std::vector<Index>          nnpAtomIndices = { 3, 4, 5, 9, 10, 11 };
    std::unique_ptr<gmx_mtop_t> mtop           = makeMtopFromFile("4water", "");

    NNPotTopologyPreprocessor topPrep(nnpAtomIndices);
    MDLogger                  logger;
    WarningHandler            wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Neural network potential Interface is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded NNP atoms: 6\nNumber of regular atoms: 6\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 6\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of settles removed: 2 \\(replaced by 4 "
                                        "InteractionFunction::ConnectBonds\\) \n");
    topPrep.preprocess(mtop.get(), logHelper_.logger(), &wi);
}

TEST_F(NNPotTopologyPreprocessorTest, AlanineDipeptideWithLinkAtomsNoConstraints)
{
    // Reference input alanine_vacuo.top
    std::vector<Index> nnpAtomIndices = { 8, 9, 10, 11, 12, 13 };
    auto               mtop           = makeMtopFromFile("alanine_vacuo", "");

    NNPotTopologyPreprocessor topPrep(nnpAtomIndices);
    MDLogger                  logger;
    WarningHandler            wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Neural network potential Interface is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded NNP atoms: 6\nNumber of regular atoms: 16\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 6\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of bonds removed: 8\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of InteractionFunction::ConnectBonds \\(type 5 bonds\\) added: 5\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of angles removed: 11\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of dihedrals removed: 9\n");
    topPrep.preprocess(mtop.get(), logHelper_.logger(), &wi);
}

TEST_F(NNPotTopologyPreprocessorTest, AlanineDipeptideWithLinkAtomsWithConstraints)
{
    // Reference input alanine_vacuo.top with constraints=all-bonds
    std::vector<Index> nnpAtomIndices = { 8, 9, 10, 11, 12, 13 };
    auto               mtop = makeMtopFromFile("alanine_vacuo", "constraints = all-bonds");

    NNPotTopologyPreprocessor topPrep(nnpAtomIndices);
    MDLogger                  logger;
    WarningHandler            wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Neural network potential Interface is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded NNP atoms: 6\nNumber of regular atoms: 16\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 6\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of bonds removed: 3\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of angles removed: 11\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of dihedrals removed: 9\n");
    topPrep.preprocess(mtop.get(), logHelper_.logger(), &wi);
}

} // namespace test

} // namespace gmx
