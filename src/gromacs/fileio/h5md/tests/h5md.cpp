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
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md.h"

#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture which sets up an empty H5md file.
using H5mdIoTest = H5mdTestBase;

/*! \brief Test that opening (creating a new), closing, re-opening and closing
 * an H5MD file works
 */
TEST(H5mdFileTest, CanCreateAndCloseH5mdFile)
{
    TestFileManager       fileManager;
    std::filesystem::path filename = fileManager.getTemporaryFilePath("ref.h5md");
    {
        EXPECT_THROW_GMX(H5md fileToRead(filename, H5mdFileMode::Read), FileIOError);
    }
    {
        gmx::H5md fileToWrite(filename, H5mdFileMode::Write);
    }
    {
        gmx::H5md fileToRead(filename, H5mdFileMode::Read);
    }
}

/*! \brief Test that writing attributes work, before closing the file and
 * after re-opening it.
 */
TEST(H5mdFileTest, CanWriteAndReadH5mdFileMetaData)
{
    TestFileManager       fileManager;
    std::filesystem::path filename = fileManager.getTemporaryFilePath("ref.h5md");
    const std::string     referenceAuthorName("AuthorName!");
    const std::string     referenceCreatorProgramName("GROMACS testing");
    const char referenceCreatorProgramVersion[] = "v. 2468"; // Testing using a char string on purpose.
    {
        SCOPED_TRACE("Testing H5MD writing.");
        gmx::H5md fileToWrite(filename, H5mdFileMode::Write);
        fileToWrite.setAuthor(referenceAuthorName);
        fileToWrite.setCreatorProgramName(referenceCreatorProgramName);
        fileToWrite.setCreatorProgramVersion(referenceCreatorProgramVersion);
        std::optional<std::string> testAuthorName = fileToWrite.author();
        ASSERT_TRUE(testAuthorName.has_value());
        EXPECT_EQ(referenceAuthorName, testAuthorName.value());
        std::optional<std::string> testCreatorProgramName = fileToWrite.creatorProgramName();
        ASSERT_TRUE(testCreatorProgramName.has_value());
        EXPECT_EQ(referenceCreatorProgramName, testCreatorProgramName.value());
        std::optional<std::string> testCreatorProgramVersion = fileToWrite.creatorProgramVersion();
        ASSERT_TRUE(testCreatorProgramVersion.has_value());
        EXPECT_EQ(referenceCreatorProgramVersion, testCreatorProgramVersion.value());
        /* It should not be possible to write an attribute that already exists. */
        EXPECT_THROW_GMX(fileToWrite.setAuthor(referenceAuthorName), FileIOError);
    }
    {
        SCOPED_TRACE("Testing H5MD reading.");
        gmx::H5md fileToRead(filename, H5mdFileMode::Read);
        {
            SCOPED_TRACE("Can't use setters on a file opened for reading");
            EXPECT_THROW_GMX(fileToRead.setAuthor(referenceAuthorName), FileIOError);
            EXPECT_THROW_GMX(fileToRead.setCreatorProgramName(referenceCreatorProgramName), FileIOError);
            EXPECT_THROW_GMX(fileToRead.setCreatorProgramVersion(referenceCreatorProgramVersion),
                             FileIOError);
        }
        std::optional<std::string> testAuthorName = fileToRead.author();
        ASSERT_TRUE(testAuthorName.has_value());
        EXPECT_EQ(referenceAuthorName, testAuthorName.value());
        std::optional<std::string> testCreatorProgramName = fileToRead.creatorProgramName();
        ASSERT_TRUE(testCreatorProgramName.has_value());
        EXPECT_EQ(referenceCreatorProgramName, testCreatorProgramName.value());
        std::optional<std::string> testCreatorProgramVersion = fileToRead.creatorProgramVersion();
        ASSERT_TRUE(testCreatorProgramVersion.has_value());
        EXPECT_EQ(referenceCreatorProgramVersion, testCreatorProgramVersion.value());
    }
}

TEST_F(H5mdIoTest, SetupFileFromInputCreatesParticlesGroup)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;

    file().setupFileFromInput(mtop, inputRecord);
    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system"));
}

TEST_F(H5mdIoTest, SetupFileFromInputThrowsForNoAtoms)
{
    gmx_mtop_t mtop;
    mtop.natoms = 0;
    t_inputrec inputRecord;
    inputRecord.nstxout = 1; // Trajectory writing is enabled for nstout >0
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 1;

    EXPECT_THROW(file().setupFileFromInput(mtop, inputRecord), gmx::FileIOError);
}

TEST_F(H5mdIoTest, SetupFileFromInputCreatesNoTrajectoryGroupsIfNoOutput)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;
    inputRecord.nstxout = 0; // Trajectory writing is not enabled for <=0
    inputRecord.nstvout = 0;
    inputRecord.nstfout = 0;

    file().setupFileFromInput(mtop, inputRecord);
    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system"));
    EXPECT_THROW(openGroup(fileid(), "/particles/system/position"), gmx::FileIOError);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/velocity"), gmx::FileIOError);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/force"), gmx::FileIOError);
}

TEST_F(H5mdIoTest, SetupFileFromInputCreatesPositionGroupIfSet)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;
    // Trajectory writing is enabled for nstout >0: only enable positions here to ensure independence
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 0;
    inputRecord.nstfout = 0;

    file().setupFileFromInput(mtop, inputRecord);
    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/position"));
    EXPECT_THROW(openGroup(fileid(), "/particles/system/velocity"), gmx::FileIOError);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/force"), gmx::FileIOError);

    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/box"))
            << "Box group must always be created";
    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/box/edges"))
            << "Edges group must only be created if position data is output";
}

TEST_F(H5mdIoTest, SetupFileFromInputCreatesVelocityGroupIfSet)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;
    // Trajectory writing is enabled for nstout >0: only enable velocities here to ensure independence
    inputRecord.nstxout = 0;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 0;

    file().setupFileFromInput(mtop, inputRecord);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/position"), gmx::FileIOError);
    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/velocity"));
    EXPECT_THROW(openGroup(fileid(), "/particles/system/force"), gmx::FileIOError);

    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/box"))
            << "Box group must always be created";
    EXPECT_THROW(openGroup(fileid(), "/particles/system/box/edges"), gmx::FileIOError)
            << "Edges group must only be created if position data is output";
}

TEST_F(H5mdIoTest, SetupFileFromInputCreatesForceGroupIfSet)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;
    // Trajectory writing is enabled for nstout >0: only enable forces here to ensure independence
    inputRecord.nstxout = 0;
    inputRecord.nstvout = 0;
    inputRecord.nstfout = 1;

    file().setupFileFromInput(mtop, inputRecord);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/position"), gmx::FileIOError);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/velocity"), gmx::FileIOError);
    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/force"));

    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/box"))
            << "Box group must always be created";
    EXPECT_THROW(openGroup(fileid(), "/particles/system/box/edges"), gmx::FileIOError)
            << "Edges group must only be created if position data is output";
}

TEST_F(H5mdIoTest, SetupFileFromInputIgnoresNstxoutCompressed)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;
    // Trajectory writing is enabled for nstout >0: only enable compressed output here
    // to assert that this does not create the position group
    inputRecord.nstxout            = 0;
    inputRecord.nstvout            = 0;
    inputRecord.nstfout            = 0;
    inputRecord.nstxout_compressed = 1;

    file().setupFileFromInput(mtop, inputRecord);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/position"), gmx::FileIOError);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/velocity"), gmx::FileIOError);
    EXPECT_THROW(openGroup(fileid(), "/particles/system/force"), gmx::FileIOError);

    EXPECT_NO_THROW(openGroup(fileid(), "/particles/system/box"))
            << "Box group must always be created";
    EXPECT_THROW(openGroup(fileid(), "/particles/system/box/edges"), gmx::FileIOError)
            << "Edges group must not be created for compressed output";
}

TEST_F(H5mdIoTest, SetupFileFromInputSetsCorrectDataSetDims)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 1;

    // Read the topology from a test system in our simulation data base
    const std::string fileNameBase = "spc2-traj";
    TprAndFileManager tprFileHandle(fileNameBase);
    bool              haveTopology;
    gmx_mtop_t        mtop;
    readConfAndTopology(tprFileHandle.tprName(), &haveTopology, &mtop, nullptr, nullptr, nullptr, nullptr);
    const hsize_t numAtoms = static_cast<hsize_t>(mtop.natoms);

    file().setupFileFromInput(mtop, inputRecord);
    const H5mdFrameDataSet<RVec> position(fileid(), "/particles/system/position/value");
    EXPECT_EQ(position.numFrames(), 0);
    EXPECT_EQ(position.frameDims(), DataSetDims{ numAtoms });
    const H5mdFrameDataSet<RVec> velocity(fileid(), "/particles/system/velocity/value");
    EXPECT_EQ(velocity.numFrames(), 0);
    EXPECT_EQ(velocity.frameDims(), DataSetDims{ numAtoms });
    const H5mdFrameDataSet<RVec> force(fileid(), "/particles/system/force/value");
    EXPECT_EQ(force.numFrames(), 0);
    EXPECT_EQ(force.frameDims(), DataSetDims{ numAtoms });
}

TEST_F(H5mdIoTest, BoxGroupForPbcXyz)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.pbcType = PbcType::Xyz;

    const hsize_t numAtoms = 6;
    gmx_mtop_t    mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);
    const auto [group, groupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/box"));

    EXPECT_EQ(getAttribute<int32_t>(group, "dimension"), DIM)
            << "Dimension attribute must be 3 for all kinds of PBC";
    EXPECT_EQ(getAttributeVector<std::string>(group, "boundary"),
              (std::vector<std::string>{ "periodic", "periodic", "periodic" }))
            << "For PBC=Xyz all boundaries are periodic";

    const H5mdFrameDataSet<real> dataSet(fileid(), "/particles/system/box/edges/value");
    EXPECT_EQ(dataSet.frameDims(), (DataSetDims{ DIM, DIM }));
}

TEST_F(H5mdIoTest, BoxGroupAttributesPbcXy)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.pbcType = PbcType::XY;

    const hsize_t numAtoms = 6;
    gmx_mtop_t    mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);
    const auto [group, groupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/box"));

    EXPECT_EQ(getAttribute<int32_t>(group, "dimension"), DIM)
            << "Dimension attribute must be 3 for all kinds of PBC";
    EXPECT_EQ(getAttributeVector<std::string>(group, "boundary"),
              std::vector<std::string>({ "periodic", "periodic", "none" }))
            << "For PBC=XY the Z value is none";

    const H5mdFrameDataSet<real> dataSet(fileid(), "/particles/system/box/edges/value");
    EXPECT_EQ(dataSet.frameDims(), (DataSetDims{ DIM, DIM }));
}

TEST_F(H5mdIoTest, BoxGroupAttributesPbcNo)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.pbcType = PbcType::No;

    const hsize_t numAtoms = 6;
    gmx_mtop_t    mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);
    const auto [group, groupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/box"));

    EXPECT_EQ(getAttribute<int32_t>(group, "dimension"), DIM)
            << "Dimension attribute must be 3 for all kinds of PBC";
    EXPECT_EQ(getAttributeVector<std::string>(group, "boundary"),
              std::vector<std::string>({ "none", "none", "none" }))
            << "For PBC=no all values are none";

    const H5mdFrameDataSet<real> dataSet(fileid(), "/particles/system/box/edges/value");
    EXPECT_EQ(dataSet.frameDims(), (DataSetDims{ DIM, DIM }));
}

TEST_F(H5mdIoTest, BoxGroupAttributesPbcScrew)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.pbcType = PbcType::Screw;

    const hsize_t numAtoms = 6;
    gmx_mtop_t    mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);
    const auto [group, groupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/box"));

    EXPECT_EQ(getAttribute<int32_t>(group, "dimension"), DIM)
            << "Dimension attribute must be 3 for all kinds of PBC";
    EXPECT_EQ(getAttributeVector<std::string>(group, "boundary"),
              std::vector<std::string>({ "none", "none", "none" }))
            << "For PBC=screw all values are none";

    const H5mdFrameDataSet<real> dataSet(fileid(), "/particles/system/box/edges/value");
    EXPECT_EQ(dataSet.frameDims(), (DataSetDims{ DIM, DIM }));
}

TEST_F(H5mdIoTest, BoxStepAndTimeDataSetsAreHardLinkedToPosition)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.pbcType = PbcType::Xyz;

    const hsize_t numAtoms = 6;
    gmx_mtop_t    mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);

    H5mdScalarFrameDataSet<int64_t> stepInPositionGroup(fileid(), "/particles/system/position/step");
    H5mdScalarFrameDataSet<double> timeInPositionGroup(fileid(), "/particles/system/position/time");

    const int64_t stepToWrite = 51;
    stepInPositionGroup.writeNextFrame(stepToWrite);
    const double timeToWrite = -501.5;
    timeInPositionGroup.writeNextFrame(timeToWrite);

    H5mdScalarFrameDataSet<int64_t> stepInBoxGroup(fileid(), "/particles/system/box/edges/step");
    int64_t                         readStepBuffer;
    stepInBoxGroup.readFrame(0, &readStepBuffer);
    EXPECT_EQ(readStepBuffer, stepToWrite)
            << "Step written to position block must exist in the box block";

    H5mdScalarFrameDataSet<double> timeInBoxGroup(fileid(), "/particles/system/box/edges/time");
    double                         readTimeBuffer;
    timeInBoxGroup.readFrame(0, &readTimeBuffer);
    EXPECT_EQ(readTimeBuffer, timeToWrite)
            << "Time written to position block must exist in the box block";
}

} // namespace
} // namespace test
} // namespace gmx
