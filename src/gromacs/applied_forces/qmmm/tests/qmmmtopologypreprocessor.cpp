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

class QMMMTopologyPreprocessorTest : public ::testing::Test
{
public:
    /*! \brief Generates tpr file from *.top and *.gro existing in the simulation database directory
     * and loads gmx_mtop_t from it
     */
    void makeMtopFromFile(const std::string& fileName, const std::string& mdpContent)
    {
        const std::filesystem::path simData =
                gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();

        // Generate empty mdp file
        const std::string mdpInputFileName =
                fileManager_.getTemporaryFilePath(fileName + ".mdp").string();
        gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpContent);

        // Generate tpr file
        const std::string tprName = fileManager_.getTemporaryFilePath(fileName + ".tpr").string();
        {
            gmx::test::CommandLine caller;
            caller.append("grompp");
            caller.addOption("-f", mdpInputFileName);
            caller.addOption("-p", (simData / fileName).replace_extension(".top").string());
            caller.addOption("-c", (simData / fileName).replace_extension(".gro").string());
            caller.addOption("-o", tprName);
            ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
        }

        // Load topology
        bool    fullTopology;
        PbcType pbcType;
        matrix  box;
        readConfAndTopology(tprName.c_str(), &fullTopology, &mtop_, &pbcType, nullptr, nullptr, box);
    }

protected:
    gmx::test::TestFileManager fileManager_;
    std::vector<Index>         qmIndices_;
    gmx_mtop_t                 mtop_;
};

class QMMMTopologyPreprocessorChecker : public gmx::test::TestReferenceChecker
{
public:
    QMMMTopologyPreprocessorChecker(const gmx::test::TestReferenceChecker& rootChecker) :
        gmx::test::TestReferenceChecker(rootChecker)
    {
        // Tolerance of all charges and vectors should be 1E-3
        setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    }

    void checkAll(const QMMMTopologyInfo& info, const gmx_mtop_t& mtop)
    {
        checkInteger(info.numMMAtoms, "Number of MM atoms");
        checkInteger(info.numQMAtoms, "Number of QM atoms");
        checkReal(info.remainingMMCharge, "MM charge");
        checkReal(info.totalClassicalChargeOfQMAtoms, "QM charge");
        checkInteger(info.numExclusionsMade, "Exclusions Made");
        checkInteger(info.numBondsRemoved, "Removed bonds");
        checkInteger(info.numAnglesRemoved, "Removed angles");
        checkInteger(info.numDihedralsRemoved, "Removed dihedrals");
        checkInteger(info.numSettleRemoved, "Removed settles");
        checkInteger(info.numConnBondsAdded, "Generated CONNBONDS");
        checkInteger(info.numVirtualSitesModified, "Removed vsites");
        checkInteger(info.numConstrainedBondsInQMSubsystem, "QM Constraints Found");
        checkInteger(info.numLinkBonds, "Link-atom sites");
        checkInteger(mtop.molblock.size(), "Molblocks in topology");
    }
};

TEST_F(QMMMTopologyPreprocessorTest, CanConstruct)
{
    qmIndices_ = { 0, 1, 2 };
    EXPECT_NO_THROW(QMMMTopologyPreprocessor topPrep(qmIndices_));
}

TEST_F(QMMMTopologyPreprocessorTest, FourWatersFirstQMNoLink)
{
    // Reference input 4x SPCE waters from database 4waters.top (first one QM) and no Link atoms
    makeMtopFromFile("4water", "");
    qmIndices_ = { 0, 1, 2 };

    QMMMTopologyPreprocessor topPrep(qmIndices_);
    topPrep.preprocess(&mtop_);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    gmx::test::TestReferenceData    data;
    QMMMTopologyPreprocessorChecker checker(data.rootChecker());
    checker.checkAll(info, mtop_);
}

TEST_F(QMMMTopologyPreprocessorTest, FourWatersSeondAndForthQMNoLink)
{
    // Reference input 4x SPCE waters from database 4waters.top (second and forth are QM) and no Link atoms
    makeMtopFromFile("4water", "");
    qmIndices_ = { 3, 4, 5, 9, 10, 11 };

    QMMMTopologyPreprocessor topPrep(qmIndices_);
    topPrep.preprocess(&mtop_);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    gmx::test::TestReferenceData    data;
    QMMMTopologyPreprocessorChecker checker(data.rootChecker());
    checker.checkAll(info, mtop_);
}

TEST_F(QMMMTopologyPreprocessorTest, FourWatersFirstQMWithLink)
{
    // Reference input 4x SPCE waters from database 4waters.top (first one QM) with Link atom
    makeMtopFromFile("4water", "");
    qmIndices_ = { 0, 1 };

    QMMMTopologyPreprocessor topPrep(qmIndices_);
    topPrep.preprocess(&mtop_);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    gmx::test::TestReferenceData    data;
    QMMMTopologyPreprocessorChecker checker(data.rootChecker());
    checker.checkAll(info, mtop_);
}

TEST_F(QMMMTopologyPreprocessorTest, AlanineDipeptideWithLinksNoConstraints)
{
    // Reference input alanine_vacuo.top
    makeMtopFromFile("alanine_vacuo", "");
    qmIndices_ = { 8, 9, 10, 11, 12, 13 };

    QMMMTopologyPreprocessor topPrep(qmIndices_);
    topPrep.preprocess(&mtop_);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    gmx::test::TestReferenceData    data;
    QMMMTopologyPreprocessorChecker checker(data.rootChecker());
    checker.checkAll(info, mtop_);
}

TEST_F(QMMMTopologyPreprocessorTest, AlanineDipeptideWithLinksWithConstraints)
{
    // Reference input alanine_vacuo.top with constraints=all-bonds
    makeMtopFromFile("alanine_vacuo", "constraints = all-bonds");
    qmIndices_ = { 8, 9, 10, 11, 12, 13 };

    QMMMTopologyPreprocessor topPrep(qmIndices_);
    topPrep.preprocess(&mtop_);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    gmx::test::TestReferenceData    data;
    QMMMTopologyPreprocessorChecker checker(data.rootChecker());
    checker.checkAll(info, mtop_);
}

TEST_F(QMMMTopologyPreprocessorTest, RemovingQMVsites)
{
    // Reference input vistes_test.top
    makeMtopFromFile("vsite_test", "");
    qmIndices_ = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    QMMMTopologyPreprocessor topPrep(qmIndices_);
    topPrep.preprocess(&mtop_);

    // Get data about changes and check it
    QMMMTopologyInfo info = topPrep.topInfo();

    gmx::test::TestReferenceData    data;
    QMMMTopologyPreprocessorChecker checker(data.rootChecker());
    checker.checkAll(info, mtop_);
}

} // namespace gmx
