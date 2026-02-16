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

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/nnpot/nnpot.h"
#include "gromacs/applied_forces/nnpot/nnpotoptions.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
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
class NNPotTopologyPreprocessorTest : public ::testing::TestWithParam<NNPotEmbedding>
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

    //! \brief Generates a default mdp values for NNPotOptions
    static KeyValueTreeObject nnpotBuildMdpValues(const std::map<std::string, std::string>& additionalValues = {})
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(std::string(NNPotModuleInfo::sc_name) + "-active",
                                              std::string("true"));
        for (const auto& [key, value] : additionalValues)
        {
            mdpValueBuilder.rootObject().addValue(std::string(NNPotModuleInfo::sc_name) + "-" + key, value);
        }
        return mdpValueBuilder.build();
    }

    //! \brief Creates an IndexGroupsAndNames object with the given atom indices
    static IndexGroupsAndNames indexGroupsAndNames(const std::vector<int>& nnpAtomIndices)
    {
        // Create an IndexGroupsAndNames object
        std::vector<IndexGroup> indexGroups;
        // "System" is the default group name for NNPot
        indexGroups.push_back({ "System", nnpAtomIndices });
        return IndexGroupsAndNames(indexGroups);
    }

    //! \brief Helper function to create an NNPotOptions object
    NNPotOptions buildDefaultOptions(const std::vector<int>& nnpAtomIndices, WarningHandler* wi)
    {
        NNPotOptions options;
        test::fillOptionsFromMdpValues(nnpotBuildMdpValues(), &options);
        options.setLogger(logHelper_.logger());
        options.setWarninp(wi);
        options.setInputGroupIndices(indexGroupsAndNames(nnpAtomIndices));
        return options;
    }

    NNPotOptions buildOptionsWithValues(const std::vector<int>&                   nnpAtomIndices,
                                        WarningHandler*                           wi,
                                        const std::map<std::string, std::string>& additionalValues)
    {
        NNPotOptions options;
        test::fillOptionsFromMdpValues(nnpotBuildMdpValues(additionalValues), &options);
        options.setLogger(logHelper_.logger());
        options.setWarninp(wi);
        options.setInputGroupIndices(indexGroupsAndNames(nnpAtomIndices));
        return options;
    }

    //! \brief Helper function to expect a log message
    void expectLogMessage(const char* msg)
    {
        logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, msg);
    }

protected:
    LoggerTestHelper logHelper_;
};

TEST_P(NNPotTopologyPreprocessorTest, FourWatersFirstInNNPRegion)
{
    const NNPotEmbedding embedding = GetParam();
    // Reference input 4x SPCE waters from database 4waters.top
    // First water is NNP input
    std::vector<int>            nnpAtomIndices = { 0, 1, 2 };
    std::unique_ptr<gmx_mtop_t> mtop           = makeMtopFromFile("4water", "");

    // Create NNPotOptions object and set required things
    WarningHandler wi(true, 0);
    NNPotOptions   options;
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        options = buildOptionsWithValues(
                nnpAtomIndices, &wi, { { "embedding", "electrostatic-model" } });
    }
    else
    {
        options = buildDefaultOptions(nnpAtomIndices, &wi);
    }

    expectLogMessage("Neural network potential interface is active, topology was modified!");
    expectLogMessage("Number of embedded NNP atoms: 3\nNumber of regular atoms: 9\n");
    expectLogMessage("Number of exclusions made: 3\n");
    expectLogMessage(
            "Number of settles removed: 1 \\(replaced by 2 "
            "InteractionFunction::ConnectBonds\\) \n");
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        expectLogMessage("Electrostatic embedding scheme is used.\n");
        expectLogMessage(
                "Total charge of the classical system \\(before modifications\\): 0.00000");
        expectLogMessage("Classical charge removed from embedded atoms: 0.00000\n");
    }
    EXPECT_NO_THROW(options.modifyTopology(mtop.get()));
}

TEST_P(NNPotTopologyPreprocessorTest, FourWatersSecondAndFourthInNNPRegion)
{
    const NNPotEmbedding embedding = GetParam();
    // Reference input 4x SPCE waters from database 4waters.top
    // second and fourth are NNP input
    std::vector<int>            nnpAtomIndices = { 3, 4, 5, 9, 10, 11 };
    std::unique_ptr<gmx_mtop_t> mtop           = makeMtopFromFile("4water", "");

    WarningHandler wi(true, 0);
    NNPotOptions   options;
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        options = buildOptionsWithValues(
                nnpAtomIndices, &wi, { { "embedding", "electrostatic-model" } });
    }
    else
    {
        options = buildDefaultOptions(nnpAtomIndices, &wi);
    }

    expectLogMessage("Neural network potential interface is active, topology was modified!");
    expectLogMessage("Number of embedded NNP atoms: 6\nNumber of regular atoms: 6\n");
    expectLogMessage("Number of exclusions made: 6\n");
    expectLogMessage(
            "Number of settles removed: 2 \\(replaced by 4 "
            "InteractionFunction::ConnectBonds\\) \n");
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        expectLogMessage("Electrostatic embedding scheme is used.\n");
        expectLogMessage(
                "Total charge of the classical system \\(before modifications\\): 0.00000");
        expectLogMessage("Classical charge removed from embedded atoms: 0.00000\n");
    }
    EXPECT_NO_THROW(options.modifyTopology(mtop.get()));
}

TEST_P(NNPotTopologyPreprocessorTest, AlanineDipeptideWithLinkAtomsNoConstraints)
{
    const NNPotEmbedding embedding = GetParam();
    // Reference input alanine_vacuo.top
    std::vector<int> nnpAtomIndices = { 8, 9, 10, 11, 12, 13 };
    auto             mtop           = makeMtopFromFile("alanine_vacuo", "");

    // Create NNPotOptions object and set required things
    WarningHandler wi(true, 0);
    NNPotOptions   options;
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        options = buildOptionsWithValues(
                nnpAtomIndices, &wi, { { "embedding", "electrostatic-model" } });
    }
    else
    {
        options = buildDefaultOptions(nnpAtomIndices, &wi);
    }

    expectLogMessage("Neural network potential interface is active, topology was modified!");
    expectLogMessage("Number of embedded NNP atoms: 6\nNumber of regular atoms: 16\n");
    expectLogMessage("Number of exclusions made: 6\n");
    expectLogMessage("Number of bonds removed: 8\n");
    expectLogMessage("Number of InteractionFunction::ConnectBonds \\(type 5 bonds\\) added: 5\n");
    expectLogMessage("Number of angles removed: 11\n");
    expectLogMessage("Number of dihedrals removed: 9\n");
    expectLogMessage("Number of link bonds added: 2\n");
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        expectLogMessage("Electrostatic embedding scheme is used.\n");
        expectLogMessage(
                "Total charge of the classical system \\(before modifications\\): 0.00000");
        expectLogMessage("Classical charge removed from embedded atoms: 0.11440\n");
    }
    EXPECT_NO_THROW(options.modifyTopology(mtop.get()));
}

TEST_P(NNPotTopologyPreprocessorTest, AlanineDipeptideWithLinkAtomsWithConstraints)
{
    const NNPotEmbedding embedding = GetParam();
    // Reference input alanine_vacuo.top with constraints=all-bonds
    std::vector<int> nnpAtomIndices = { 8, 9, 10, 11, 12, 13 };
    auto             mtop           = makeMtopFromFile("alanine_vacuo", "constraints = all-bonds");

    // Create NNPotOptions object and set required things
    WarningHandler wi(true, 0);
    NNPotOptions   options;
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        options = buildOptionsWithValues(
                nnpAtomIndices, &wi, { { "embedding", "electrostatic-model" } });
    }
    else
    {
        options = buildDefaultOptions(nnpAtomIndices, &wi);
    }

    expectLogMessage("Neural network potential interface is active, topology was modified!");
    expectLogMessage("Number of embedded NNP atoms: 6\nNumber of regular atoms: 16\n");
    expectLogMessage("Number of exclusions made: 6\n");
    expectLogMessage("Number of bonds removed: 3\n");
    expectLogMessage("Number of angles removed: 11\n");
    expectLogMessage("Number of dihedrals removed: 9\n");
    expectLogMessage("Number of link bonds added: 2\n");
    if (GetParam() == NNPotEmbedding::ElectrostaticModel)
    {
        expectLogMessage("Electrostatic embedding scheme is used.\n");
        expectLogMessage(
                "Total charge of the classical system \\(before modifications\\): 0.00000");
        expectLogMessage("Classical charge removed from embedded atoms: 0.11440\n");
    }
    EXPECT_NO_THROW(options.modifyTopology(mtop.get()));

    // expect one warning about constrained bonds
    int numWarnings = embedding == NNPotEmbedding::ElectrostaticModel ? 2 : 1;
    ASSERT_EQ(wi.warningCount(), numWarnings);
}

TEST_P(NNPotTopologyPreprocessorTest, RemovingChargeOnVSites)
{
    const NNPotEmbedding embedding = GetParam();
    // Reference input vistes_test.top
    std::unique_ptr<gmx_mtop_t> mtop = makeMtopFromFile("vsite_test", "");
    std::vector<int> qmIndices       = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    WarningHandler wi(true, 0);
    NNPotOptions   options;
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        options = buildOptionsWithValues(qmIndices, &wi, { { "embedding", "electrostatic-model" } });
    }
    else
    {
        options = buildDefaultOptions(qmIndices, &wi);
    }

    expectLogMessage("Neural network potential interface is active, topology was modified!");
    expectLogMessage("Number of embedded NNP atoms: 16\nNumber of regular atoms: 8\n");
    expectLogMessage("Number of exclusions made: 16\n");
    expectLogMessage("Number of bonds removed: 28\n");
    expectLogMessage("Number of InteractionFunction::ConnectBonds \\(type 5 bonds\\) added: 15\n");
    expectLogMessage("Number of angles removed: 14\n");
    expectLogMessage("Number of dihedrals removed: 13\n");
    if (embedding == NNPotEmbedding::ElectrostaticModel)
    {
        expectLogMessage("Electrostatic embedding scheme is used.\n");
        expectLogMessage(
                "Total charge of the classical system \\(before modifications\\): 0.00000");
        expectLogMessage("Classical charge removed from embedded atoms: 0.00000\n");
        expectLogMessage(
                "Note: There are 8 virtual sites found, which are built from embedded atoms only. "
                "Classical charges on them have been removed as well.\n");
    }
    EXPECT_NO_THROW(options.modifyTopology(mtop.get()));
}

INSTANTIATE_TEST_SUITE_P(NNPotTopologyPreprocessorTests,
                         NNPotTopologyPreprocessorTest,
                         ::testing::Values(NNPotEmbedding::Mechanical, NNPotEmbedding::ElectrostaticModel));

} // namespace test

} // namespace gmx
