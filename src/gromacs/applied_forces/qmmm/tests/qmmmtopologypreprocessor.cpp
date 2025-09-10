/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for QMMMInputGenerator class for QMMM MDModule
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/qmmm/qmmmtopologypreprocessor.h"

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

class QMMMTopologyPreprocessorTest : public ::testing::Test
{
public:
    /*! \brief Generates tpr file from *.top and *.gro existing in the simulation database directory
     * and loads gmx_mtop_t from it
     */
    static std::unique_ptr<gmx_mtop_t> makeMtopFromFile(const std::string& fileName,
                                                        const std::string& mdpContent)
    {
        const std::filesystem::path simData =
                gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();
        TestFileManager fileManager;

        // Generate empty mdp file
        const std::string mdpInputFileName = fileManager.getTemporaryFilePath(fileName + ".mdp").string();
        gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpContent);

        // Generate tpr file
        const std::string tprName = fileManager.getTemporaryFilePath(fileName + ".tpr").string();
        {
            gmx::test::CommandLine caller;
            caller.append("grompp");
            caller.addOption("-f", mdpInputFileName);
            caller.addOption("-p", (simData / fileName).replace_extension(".top").string());
            caller.addOption("-c", (simData / fileName).replace_extension(".gro").string());
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

TEST_F(QMMMTopologyPreprocessorTest, CanConstruct)
{
    std::vector<Index> qmIndices = { 0, 1, 2 };
    EXPECT_NO_THROW(QMMMTopologyPreprocessor topPrep(qmIndices));
}

TEST_F(QMMMTopologyPreprocessorTest, FourWatersFirstQMNoLink)
{
    // Reference input 4x SPCE waters from database 4waters.top (first one QM) and no Link atoms
    std::unique_ptr<gmx_mtop_t> mtop      = makeMtopFromFile("4water", "");
    std::vector<Index>          qmIndices = { 0, 1, 2 };

    QMMMTopologyPreprocessor topPrep(qmIndices);
    WarningHandler           wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info, "QMMM Interface with CP2K is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded QM atoms: 3\nNumber of regular atoms: 9\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Total charge of the classical system \\(before modifications\\): 0.00000");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Classical charge removed from embedded atoms: 0.00000\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 3\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of settles removed: 1 \\(replaced by 2 F_CONNBONDS\\) \n");
    topPrep.preprocess(mtop.get(), 0.0, logHelper_.logger(), &wi);
}

TEST_F(QMMMTopologyPreprocessorTest, FourWatersSeondAndForthQMNoLink)
{
    // Reference input 4x SPCE waters from database 4waters.top (second and forth are QM) and no Link atoms
    std::unique_ptr<gmx_mtop_t> mtop      = makeMtopFromFile("4water", "");
    std::vector<Index>          qmIndices = { 3, 4, 5, 9, 10, 11 };

    QMMMTopologyPreprocessor topPrep(qmIndices);
    WarningHandler           wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info, "QMMM Interface with CP2K is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded QM atoms: 6\nNumber of regular atoms: 6\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Total charge of the classical system \\(before modifications\\): 0.00000");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Classical charge removed from embedded atoms: 0.00000\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 6\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of settles removed: 2 \\(replaced by 4 F_CONNBONDS\\) \n");
    topPrep.preprocess(mtop.get(), 0.0, logHelper_.logger(), &wi);
}

TEST_F(QMMMTopologyPreprocessorTest, FourWatersFirstQMWithLink)
{
    // Reference input 4x SPCE waters from database 4waters.top (first one QM) with Link atom
    std::unique_ptr<gmx_mtop_t> mtop      = makeMtopFromFile("4water", "");
    std::vector<Index>          qmIndices = { 0, 1 };

    QMMMTopologyPreprocessor topPrep(qmIndices);
    MDLogger                 logger;
    WarningHandler           wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info, "QMMM Interface with CP2K is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded QM atoms: 2\nNumber of regular atoms: 10\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Total charge of the classical system \\(before modifications\\): 0.00000");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Classical charge removed from embedded atoms: -0.41000\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 2\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of settles removed: 1 \\(replaced by 2 F_CONNBONDS\\) \n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of link bonds added: 1\n");
    topPrep.preprocess(mtop.get(), -0.41, logHelper_.logger(), &wi);
}

TEST_F(QMMMTopologyPreprocessorTest, AlanineDipeptideWithLinksNoConstraints)
{
    // Reference input alanine_vacuo.top
    std::unique_ptr<gmx_mtop_t> mtop      = makeMtopFromFile("alanine_vacuo", "");
    std::vector<Index>          qmIndices = { 8, 9, 10, 11, 12, 13 };

    QMMMTopologyPreprocessor topPrep(qmIndices);
    MDLogger                 logger;
    WarningHandler           wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info, "QMMM Interface with CP2K is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded QM atoms: 6\nNumber of regular atoms: 16\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Total charge of the classical system \\(before modifications\\): 0.00000");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Classical charge removed from embedded atoms: 0.11440\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 6\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of bonds removed: 8\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of F_CONNBONDS \\(type 5 bonds\\) added: 5\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of angles removed: 11\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of dihedrals removed: 9\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of link bonds added: 2\n");
    topPrep.preprocess(mtop.get(), 0.1144, logHelper_.logger(), &wi);
}

TEST_F(QMMMTopologyPreprocessorTest, AlanineDipeptideWithLinksWithConstraints)
{
    // Reference input alanine_vacuo.top with constraints=all-bonds
    std::unique_ptr<gmx_mtop_t> mtop = makeMtopFromFile("alanine_vacuo", "constraints = all-bonds");
    std::vector<Index>          qmIndices = { 8, 9, 10, 11, 12, 13 };

    QMMMTopologyPreprocessor topPrep(qmIndices);
    MDLogger                 logger;
    WarningHandler           wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info, "QMMM Interface with CP2K is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded QM atoms: 6\nNumber of regular atoms: 16\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Total charge of the classical system \\(before modifications\\): 0.00000");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Classical charge removed from embedded atoms: 0.11440\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 6\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of bonds removed: 3\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of angles removed: 11\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of dihedrals removed: 9\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of link bonds added: 2\n");
    topPrep.preprocess(mtop.get(), 0.1144, logHelper_.logger(), &wi);
}

TEST_F(QMMMTopologyPreprocessorTest, RemovingQMVsites)
{
    // Reference input vistes_test.top
    std::unique_ptr<gmx_mtop_t> mtop = makeMtopFromFile("vsite_test", "");
    std::vector<Index> qmIndices     = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    QMMMTopologyPreprocessor topPrep(qmIndices);
    MDLogger                 logger;
    WarningHandler           wi(true, 0);
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info, "QMMM Interface with CP2K is active, topology was modified!");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Number of embedded QM atoms: 16\nNumber of regular atoms: 8\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Total charge of the classical system \\(before modifications\\): 0.00000");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Classical charge removed from embedded atoms: 0.00000\n");
    logHelper_.expectEntryMatchingRegex(
            MDLogger::LogLevel::Info,
            "Note: There are 8 virtual sites found, which are built from embedded atoms only. "
            "Classical charges on them have been removed as well.\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of exclusions made: 16\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of bonds removed: 28\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of F_CONNBONDS \\(type 5 bonds\\) added: 15\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info, "Number of angles removed: 14\n");
    logHelper_.expectEntryMatchingRegex(MDLogger::LogLevel::Info,
                                        "Number of dihedrals removed: 13\n");
    topPrep.preprocess(mtop.get(), 0.0, logHelper_.logger(), &wi);
}

} // namespace test

} // namespace gmx
