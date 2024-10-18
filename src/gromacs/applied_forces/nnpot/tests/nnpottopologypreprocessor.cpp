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
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
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
     * \param[out] mtop          gmx_mtop_t struct to load the topology into
     */
    static void makeMtopFromFile(const std::string& simulationName,
                                 const std::string& mdpContent,
                                 gmx_mtop_t*        mtop)
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
            ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
        }

        // Load topology
        bool    fullTopology;
        PbcType pbcType;
        matrix  box;
        readConfAndTopology(tprName.c_str(), &fullTopology, mtop, &pbcType, nullptr, nullptr, box);
    }

    /*! \brief Helper function to check the processed topology.
     *
     * \param[in] info    QMMMTopologyInfo struct to check
     * \param[in] mtop    gmx_mtop_t struct to check modified molblocks
     * \param[in] checker TestReferenceChecker to use for checking
     */
    static void checkTopologyInfo(const QMMMTopologyInfo& info,
                                  const gmx_mtop_t&       mtop,
                                  TestReferenceChecker*   checker)
    {
        // Tolerance of all charges and vectors should be 1E-3
        checker->setDefaultTolerance(absoluteTolerance(0.001));

        checker->checkInteger(info.numMMAtoms, "Number of MM atoms");
        checker->checkInteger(info.numQMAtoms, "Number of NNP atoms");
        checker->checkReal(info.remainingMMCharge, "MM charge");
        checker->checkReal(info.totalClassicalChargeOfQMAtoms, "NNP charge");
        checker->checkInteger(info.numExclusionsMade, "Exclusions Made");
        checker->checkInteger(info.numBondsRemoved, "Removed bonds");
        checker->checkInteger(info.numAnglesRemoved, "Removed angles");
        checker->checkInteger(info.numDihedralsRemoved, "Removed dihedrals");
        checker->checkInteger(info.numSettleRemoved, "Removed settles");
        checker->checkInteger(info.numConnBondsAdded, "Generated CONNBONDS");
        checker->checkInteger(info.numVirtualSitesModified, "Removed vsites");
        checker->checkInteger(mtop.molblock.size(), "Molblocks in topology");
    }
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
    std::vector<Index> nnpAtomIndices = { 0, 1, 2 };
    gmx_mtop_t         mtop;
    makeMtopFromFile("4water", "", &mtop);

    NNPotTopologyPreprocessor topPrep(nnpAtomIndices);
    topPrep.preprocess(&mtop);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    gmx::test::TestReferenceData data;
    TestReferenceChecker         checker(data.rootChecker());
    checkTopologyInfo(info, mtop, &checker);
}

TEST_F(NNPotTopologyPreprocessorTest, FourWatersSecondAndFourthInQMRegion)
{
    // Reference input 4x SPCE waters from database 4waters.top
    // second and fourth are NNP input
    std::vector<Index> nnpAtomIndices = { 3, 4, 5, 9, 10, 11 };
    gmx_mtop_t         mtop;
    makeMtopFromFile("4water", "", &mtop);

    NNPotTopologyPreprocessor topPrep(nnpAtomIndices);
    topPrep.preprocess(&mtop);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checkTopologyInfo(info, mtop, &checker);
}

TEST_F(NNPotTopologyPreprocessorTest, AlanineDipeptideWithLinkAtomsNoConstraints)
{
    // Reference input alanine_vacuo.top
    std::vector<Index> nnpAtomIndices = { 8, 9, 10, 11, 12, 13 };
    gmx_mtop_t         mtop;
    makeMtopFromFile("alanine_vacuo", "", &mtop);

    QMMMTopologyPreprocessor topPrep(nnpAtomIndices);
    topPrep.preprocess(&mtop);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checkTopologyInfo(info, mtop, &checker);
}

TEST_F(NNPotTopologyPreprocessorTest, AlanineDipeptideWithLinkAtomsWithConstraints)
{
    // Reference input alanine_vacuo.top with constraints=all-bonds
    std::vector<Index> nnpAtomIndices = { 8, 9, 10, 11, 12, 13 };
    gmx_mtop_t         mtop;
    makeMtopFromFile("alanine_vacuo", "constraints = all-bonds", &mtop);

    NNPotTopologyPreprocessor topPrep(nnpAtomIndices);
    topPrep.preprocess(&mtop);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checkTopologyInfo(info, mtop, &checker);
}

} // namespace test

} // namespace gmx
