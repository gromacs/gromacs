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
 * Tests for H5MD file I/O routines
 *
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/h5md_io.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"


// #include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"


namespace
{
/*! \brief
 * Convienience type for testing combinations of atoms and number of frames.
 *
 * Fields are: atomCount, numFrames, xCompressionPrecision
 */
using H5mdWriteReadTrajectoryTestParams = std::tuple<int, int, float>;

class H5mdIoTest : public ::testing::Test, public ::testing::WithParamInterface<H5mdWriteReadTrajectoryTestParams>
{
public:
    H5mdIoTest()
    {
        clear_mat(refBox_);
        referenceFilename_ = fileManager_.getTemporaryFilePath(getFileSuffix("ref")).u8string();
        open_symtab(&topologySymbolTable_);
        strcpy(atomNameBase, "мệ🚀");
        refX_ = nullptr;
        refV_ = nullptr;
        refF_ = nullptr;
    }

    ~H5mdIoTest() override
    {
        done_symtab(&topologySymbolTable_);
        sfree(refX_);
        sfree(refV_);
        sfree(refF_);
    }

    /*! \brief Open a file used as reference for further tests. */
    void openReferenceFile(const char mode) { referenceH5mdIo_.openFile(referenceFilename_, mode); }

    /*! \brief Close the reference file. */
    void closeReferenceFile() { referenceH5mdIo_.closeFile(); }

    /*! \brief Check whether the reference file is open.
     * \returns true if the reference file is open, otherwise false. */
    bool isReferenceFileOpen() { return referenceH5mdIo_.isFileOpen(); }

    /*! \brief Set the number of reference atoms to use. */
    void setRefAtomCount(int atomCount) { refAtomCount_ = atomCount; }

    /*! \brief Set the number of reference atoms. */
    int getRefAtomCount() { return refAtomCount_; }

    void checkH5mdRootVersionNumber()
    {
        std::string fileVersion   = referenceH5mdIo_.getH5mdRootVersionNumber();
        std::string versionString = std::to_string(gmx::h5mdio::c_h5mdMajorVersion) + "."
                                    + std::to_string(gmx::h5mdio::c_h5mdMinorVersion);
        EXPECT_STREQ(fileVersion.c_str(), versionString.c_str());
    }

    void setAuthorAndCreator(std::string authorName, std::string creatorProgramName, std::string creatorProgramVersion)
    {
        referenceH5mdIo_.setAuthor(authorName);
        referenceH5mdIo_.setCreatorProgramName(creatorProgramName);
        referenceH5mdIo_.setCreatorProgramVersion(creatorProgramVersion);
    }

    void checkAuthorAndCreator(std::string referenceAuthorName,
                               std::string referenceCreatorProgramName,
                               std::string referenceCreatorProgramVersion)
    {
        std::string authorName = referenceH5mdIo_.getAuthor();
        EXPECT_STREQ(referenceAuthorName.c_str(), authorName.c_str());
        std::string creatorProgramName = referenceH5mdIo_.getCreatorProgramName();
        EXPECT_STREQ(referenceCreatorProgramName.c_str(), creatorProgramName.c_str());
        std::string creatorProgramVersion = referenceH5mdIo_.getCreatorProgramVersion();
        EXPECT_STREQ(referenceCreatorProgramVersion.c_str(), creatorProgramVersion.c_str());
    }

    /*! \brief Set the lossy compression precision to use when writing reference data. */
    void setRefCompressionPrecision(real compressionPrecision)
    {
        refCompressionPrecision_ = compressionPrecision;
    }

    /*! \brief Get the lossy compression precision of the referece data. */
    real getRefCompressionPrecision() { return refCompressionPrecision_; }

    /*! \brief Generate reference data (coordinates and velocities) for the reference atoms. */
    void generateReferenceCoordinatesAndVelocities(int frame = 0)
    {
        clear_mat(refBox_);
        refBox_[XX][XX] = 2;
        refBox_[YY][YY] = 3;
        refBox_[ZZ][ZZ] = 4.123;
        sfree(refX_);
        sfree(refV_);
        sfree(refF_);
        snew(refX_, refAtomCount_);
        snew(refV_, refAtomCount_);
        refF_       = nullptr;
        int divisor = refAtomCount_ % 8 == 0 ? 4 : 3;
        for (size_t i = 0; i < refAtomCount_; ++i)
        {
            gmx::RVec v(0.1, 0.22, -frame * 0.001);
            gmx::RVec x(i / divisor + 0.05 * (i % divisor) + frame * v[0],
                        i / divisor + 0.05 * (i % divisor) + frame * v[1],
                        i / divisor + 0.05 * (i % divisor) + frame * v[2]);
            copy_rvec(x, refX_[i]);
            copy_rvec(v, refV_[i]);
        }
    }

    /* Initialize the molecular system information in the reference H5MD file. */
    void setupMolecularSystem()
    {
        gmx::ArrayRef<const int> index;
        std::string              indexGroupName = "";
        gmx::setupMolecularSystem(&referenceH5mdIo_, refTopology_, index, indexGroupName);
    }

    /* Generate a bunch of nonsensical atom names for the reference data, using UTF8 characters. */
    void generateReferenceTopologyAtomNames()
    {
        auto& moltype    = refTopology_.moltype.emplace_back();
        moltype.atoms.nr = refAtomCount_;
        init_t_atoms(&moltype.atoms, refAtomCount_, false);

        /* We are using some UTF8 characters, so there must be some margin in the string length. */
        for (size_t i = 0; i < refAtomCount_; i++)
        {
            char tmpUtf8CharBuffer[gmx::h5mdio::c_atomStringLen + 1];
            sprintf(tmpUtf8CharBuffer, "%s%zu", atomNameBase, i % 1000);
            moltype.atoms.atom[i].resind = 0;
            moltype.atoms.atomname[i]    = put_symtab(&topologySymbolTable_, tmpUtf8CharBuffer);
        }

        refTopology_.molblock.resize(1);
        refTopology_.molblock[0].type = 0;
        refTopology_.molblock[0].nmol = 1;
        refTopology_.natoms           = moltype.atoms.nr * refTopology_.molblock[0].nmol;

        refTopology_.finalize();
    }

    std::vector<std::string> readAtomNamesFromReferenceFile()
    {
        return referenceH5mdIo_.readStringProperty("/particles/system", "atomname");
    }

    void compareAtomNamesToReference(const std::vector<std::string>& atomNames)
    {
        EXPECT_EQ(refAtomCount_, atomNames.size());
        for (size_t i = 0; i < refAtomCount_; i++)
        {
            char tmpUtf8CharBuffer[gmx::h5mdio::c_atomStringLen + 1];
            sprintf(tmpUtf8CharBuffer, "%s%zu", atomNameBase, i % 1000);
            // printf("name %zu: %s %s\n", i, atomNames[i].c_str(), tmpUtf8CharBuffer);
            EXPECT_STREQ(tmpUtf8CharBuffer, atomNames[i].c_str());
        }
    }

    void generateReferenceTopology() {}

    void writeReferenceTrajectoryFrame(int step, real time, real lambda)
    {
        gmx::writeFrameToStandardDataBlocks(
                &referenceH5mdIo_, step, time, lambda, refBox_, refAtomCount_, refX_, refV_, refF_, refCompressionPrecision_);
    }

    int64_t readReferenceNumAtoms(const std::string dataBlockName)
    {
        return referenceH5mdIo_.getNumberOfParticles(dataBlockName);
    }

    int64_t readReferenceNumFrames(const std::string dataBlockName)
    {
        return referenceH5mdIo_.getNumberOfFrames(dataBlockName);
    }

    real getFirstTimeFromAllDataBlocks()
    {
        return referenceH5mdIo_.getFirstTimeFromAllDataBlocks();
    }

    real getFinalTimeFromAllDataBlocks()
    {
        return referenceH5mdIo_.getFinalTimeFromAllDataBlocks();
    }

    void readNextFrameAndCompareToReference(int64_t referenceStep, int64_t referenceNumFrames)
    {
        int64_t testStep;
        real    testTime;
        real    testLambda;
        matrix  testBox;
        real    testPrecision;
        bool    testReadLambda, testReadBox, testReadX, testReadV, testReadF;
        rvec*   testX;
        rvec*   testV;
        rvec*   testF;
        snew(testX, refAtomCount_);
        snew(testV, refAtomCount_);
        snew(testF, refAtomCount_);
        gmx::readNextFrameOfStandardDataBlocks(&referenceH5mdIo_,
                                               &testStep,
                                               &testTime,
                                               &testLambda,
                                               testBox,
                                               testX,
                                               testV,
                                               testF,
                                               &testPrecision,
                                               &testReadLambda,
                                               &testReadBox,
                                               &testReadX,
                                               &testReadV,
                                               &testReadF);

        real referenceTime   = referenceStep * 10;
        real referenceLambda = static_cast<real>(referenceStep) / referenceNumFrames;

        EXPECT_EQ(referenceStep, testStep);
        EXPECT_REAL_EQ_TOL(referenceTime, testTime, gmx::test::defaultRealTolerance());
        EXPECT_TRUE(testReadLambda);
        EXPECT_TRUE(testReadBox);
        EXPECT_TRUE(testReadX);
        EXPECT_TRUE(testReadV);
        EXPECT_FALSE(testReadF);
        EXPECT_REAL_EQ_TOL(referenceLambda, testLambda, gmx::test::defaultRealTolerance());

        if (refCompressionPrecision_ == 0)
        {

            EXPECT_EQ(-1, testPrecision);
        }
        else
        {
            EXPECT_REAL_EQ_TOL(refCompressionPrecision_, testPrecision, gmx::test::defaultRealTolerance());
        }
        for (int d1 = 0; d1 < DIM; d1++)
        {
            for (int d2 = 0; d2 < DIM; d2++)
            {
                EXPECT_REAL_EQ_TOL(refBox_[d1][d2], testBox[d1][d2], gmx::test::defaultRealTolerance());
            }
        }
        for (size_t atom = 0; atom < refAtomCount_; atom++)
        {
            for (int d = 0; d < DIM; d++)
            {
                if (refCompressionPrecision_ > 0)
                {
                    EXPECT_NEAR(refX_[atom][d], testX[atom][d], refCompressionPrecision_);
                }
                else
                {
                    EXPECT_REAL_EQ_TOL(refX_[atom][d], testX[atom][d], gmx::test::defaultRealTolerance());
                }
                EXPECT_REAL_EQ_TOL(refV_[atom][d], testV[atom][d], gmx::test::defaultRealTolerance());
            }
        }
        sfree(testX);
        sfree(testV);
        sfree(testF);
    }

    void setupWaterAndLigandGroups()
    {
        refWaterAtomNames_       = { "OW", "HW1", "HW2", "OW", "HW1", "HW2", "OW", "HW1", "HW2" };
        refWaterPartialCharges_  = { -0.83,  -0.415, -0.415, -0.83, -0.415,
                                    -0.415, -0.83,  -0.415, -0.415 };
        refLigandAtomNames_      = { "C1", "HC11", "HC12", "HC13", "C2", "HC21", "HC22", "HC23" };
        refLigandPartialCharges_ = { -0.27, 0.09, 0.09, 0.09, -0.27, 0.09, 0.09, 0.09 };

        GMX_ASSERT(refWaterAtomNames_.size() == refWaterPartialCharges_.size()
                           && refLigandAtomNames_.size() == refLigandPartialCharges_.size(),
                   "Mismatching number of elements");
    }

    void writeWaterAndLigandGroups()
    {
        /* Test variable-length string writing. */
        referenceH5mdIo_.setStringProperty("/particles/water", "atomname", refWaterAtomNames_, false, 0);
        referenceH5mdIo_.setStringProperty("/particles/ligand", "atomname", refLigandAtomNames_, false);
        referenceH5mdIo_.setFloatProperty("/particles/water", "charge", refWaterPartialCharges_, false);
        referenceH5mdIo_.setFloatProperty("/particles/ligand", "charge", refLigandPartialCharges_, false);
    }

    void readAndCheckWaterAndLigandGroups()
    {
        std::vector<std::string> testWaterAtomNames =
                referenceH5mdIo_.readStringProperty("/particles/water", "atomname");
        std::vector<std::string> testLigandAtomNames =
                referenceH5mdIo_.readStringProperty("/particles/ligand", "atomname");
        std::vector<real> testWaterCharges =
                referenceH5mdIo_.readFloatProperty("/particles/water", "charge");
        std::vector<real> testLigandCharges =
                referenceH5mdIo_.readFloatProperty("/particles/ligand", "charge");

        size_t refNumWaterNameElements  = refWaterAtomNames_.size();
        size_t refNumLigandNameElements = refLigandAtomNames_.size();

        EXPECT_EQ(refNumWaterNameElements, testWaterAtomNames.size());
        EXPECT_EQ(refNumWaterNameElements, testWaterCharges.size());
        EXPECT_EQ(refNumLigandNameElements, testLigandAtomNames.size());
        EXPECT_EQ(refNumLigandNameElements, testLigandCharges.size());

        for (size_t i = 0; i < refNumWaterNameElements; i++)
        {
            EXPECT_STREQ(refWaterAtomNames_[i].c_str(), testWaterAtomNames[i].c_str());
            EXPECT_REAL_EQ_TOL(
                    refWaterPartialCharges_[i], testWaterCharges[i], gmx::test::defaultRealTolerance());
        }
        for (size_t i = 0; i < refNumLigandNameElements; i++)
        {
            EXPECT_STREQ(refLigandAtomNames_[i].c_str(), testLigandAtomNames[i].c_str());
            EXPECT_REAL_EQ_TOL(
                    refLigandPartialCharges_[i], testLigandCharges[i], gmx::test::defaultRealTolerance());
        }
    }

private:
    static std::string getFileSuffix(const char* type)
    {
        return std::string(type) + "." + ftp2ext(efH5MD);
    }


    gmx::test::TestFileManager fileManager_;
    std::string                referenceFilename_;
    gmx::h5mdio::GmxH5mdIo     referenceH5mdIo_;
    rvec*                      refX_;
    rvec*                      refV_;
    rvec*                      refF_;
    matrix                     refBox_;
    gmx_mtop_t                 refTopology_;
    t_symtab                   topologySymbolTable_;
    size_t                     refAtomCount_;
    real                       refCompressionPrecision_;
    char                       atomNameBase[gmx::h5mdio::c_atomStringLen - 5];
    std::vector<std::string>   refWaterAtomNames_;
    std::vector<real>          refWaterPartialCharges_;
    std::vector<std::string>   refLigandAtomNames_;
    std::vector<real>          refLigandPartialCharges_;
};

/*! \brief Tests that opening (creating a new), closing, re-opening and closing
 * an H5MD file works
 */
TEST_F(H5mdIoTest, CanCreateAndCloseH5mdFile)
{
    EXPECT_THROW_GMX(openReferenceFile('r'), gmx::FileIOError);
    EXPECT_FALSE(isReferenceFileOpen());
    openReferenceFile('w');
    EXPECT_TRUE(isReferenceFileOpen());
    const std::string referenceAuthorName("AuthorName!");
    const std::string referenceCreatorProgramName("GROMACS testing");
    const char referenceCreatorProgramVersion[] = "v. 2468"; // Testing using a char string on purpose.
    setAuthorAndCreator(referenceAuthorName, referenceCreatorProgramName, referenceCreatorProgramVersion);
    checkAuthorAndCreator(referenceAuthorName, referenceCreatorProgramName, referenceCreatorProgramVersion);
    checkH5mdRootVersionNumber();
    closeReferenceFile();
    EXPECT_FALSE(isReferenceFileOpen());
    openReferenceFile('r');
    /* Verify that the version number is the same after flushing, closing and opening. */
    checkH5mdRootVersionNumber();
    EXPECT_TRUE(isReferenceFileOpen());
    closeReferenceFile();
    EXPECT_FALSE(isReferenceFileOpen());
}

TEST_F(H5mdIoTest, SelectionGroups)
{
    openReferenceFile('w');
    setupWaterAndLigandGroups();
    writeWaterAndLigandGroups();
    closeReferenceFile();
    openReferenceFile('r');
    readAndCheckWaterAndLigandGroups();
    closeReferenceFile();
}

TEST_P(H5mdIoTest, HighLevelWriteRead)
{
    auto params = GetParam();
    setRefAtomCount(std::get<0>(params));
    int numFrames = std::get<1>(params);
    setRefCompressionPrecision(std::get<2>(params));

    generateReferenceTopologyAtomNames();

    EXPECT_FALSE(isReferenceFileOpen());
    openReferenceFile('w');
    setupMolecularSystem();

    for (int i = 0; i < numFrames; i++)
    {
        real time   = i * 10;
        real lambda = static_cast<real>(i) / numFrames;
        generateReferenceCoordinatesAndVelocities(i);
        if (getRefAtomCount() <= 0)
        {
            EXPECT_THROW_GMX(writeReferenceTrajectoryFrame(i, time, lambda), gmx::FileIOError);
        }
        else
        {
            writeReferenceTrajectoryFrame(i, time, lambda);
        }
    }
    closeReferenceFile();

    openReferenceFile('r');
    EXPECT_TRUE(isReferenceFileOpen());
    if (getRefAtomCount() <= 0)
    {
        /* The rest of the tests will fail. */
        return;
    }
    std::vector<std::string> atomNames = readAtomNamesFromReferenceFile();
    compareAtomNamesToReference(atomNames);

    EXPECT_EQ(getRefAtomCount(), readReferenceNumAtoms("position"));
    EXPECT_EQ(getRefAtomCount(), readReferenceNumAtoms("velocity"));
    EXPECT_EQ(-1, readReferenceNumAtoms("force"));
    EXPECT_EQ(numFrames, readReferenceNumFrames("position"));
    EXPECT_EQ(numFrames, readReferenceNumFrames("velocity"));
    // EXPECT_THROW_GMX(readReferenceNumFrames("force"), gmx::FileIOError);
    EXPECT_EQ(0, getFirstTimeFromAllDataBlocks());
    EXPECT_EQ((numFrames - 1) * 10, getFinalTimeFromAllDataBlocks());

    for (int i = 0; i < numFrames; i++)
    {
        generateReferenceCoordinatesAndVelocities(i);

        readNextFrameAndCompareToReference(i, numFrames);
    }
    closeReferenceFile();
    EXPECT_FALSE(isReferenceFileOpen());
}

INSTANTIATE_TEST_SUITE_P(H5mdTestWriteReadCombinations,
                         H5mdIoTest,
                         ::testing::Combine(::testing::Values(0, 2, 18, 48),
                                            ::testing::Values(2, 20, 101),
                                            ::testing::Values(0, 0.001, 0.0123))); // FIXME: Why does 0.0123 cause segfaults


} // namespace
