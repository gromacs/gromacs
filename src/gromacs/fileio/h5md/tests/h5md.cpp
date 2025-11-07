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
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture which sets up an empty H5md file.
using H5mdIoTest = H5mdTestBase;

//! \brief Dummy simulation box, used in tests where this is not needed.
constexpr matrix c_unusedBox = { { -1.0, -1.0, -1.0 }, { -1.0, -1.0, -1.0 }, { -1.0, -1.0, -1.0 } };

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

TEST_F(H5mdIoTest, SetupFileFromInputWritesMetadataGroup)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;

    file().setupFileFromInput(mtop, inputRecord);

    const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(fileid(), "h5md"));

    EXPECT_TRUE(getAttributeVector<int>(group, "version").has_value())
            << "H5md specification version attribute must be written";
    EXPECT_EQ(getAttributeVector<int>(group, "version").value(), (std::vector<int>{ 1, 1 }))
            << "H5md specification version should be 1.1";

    const auto [authorGroup, authorGroupGuard] = makeH5mdGroupGuard(openGroup(group, "author"));
    EXPECT_TRUE(getAttribute<std::string>(authorGroup, "name").has_value())
            << "Author name must be written";

    const auto [creatorGroup, creatorGroupGuard] = makeH5mdGroupGuard(openGroup(group, "creator"));
    EXPECT_EQ(getAttribute<std::string>(creatorGroup, "name").value_or(""), "GROMACS")
            << "Creator name must be GROMACS";

    EXPECT_EQ(getAttribute<std::string>(creatorGroup, "version").value_or(""), gmx_version())
            << "Program version must match";
}

TEST_F(H5mdIoTest, SetupFileFromInputWritesModuleInformation)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;
    file().setupFileFromInput(mtop, inputRecord);

    const auto [metadataGroup, metadataGroupGuard] = makeH5mdGroupGuard(openGroup(fileid(), "h5md"));
    ASSERT_TRUE(objectExists(metadataGroup, "modules"))
            << "modules group must exist in /h5md metadata group after setup";

    // Check that our GROMACS module exists and has an experimental version
    const auto [gromacsGroup, gromacsGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/h5md/modules/gromacs"));
    EXPECT_EQ(getAttributeVector<int>(gromacsGroup, "version").value()[0], 0)
            << "GROMACS module specification must have experimental version 0.x to follow semantic "
               "versioning";
    // Ideally we would test some of its content here, but until that is specified
    // we only check its existence

    // Check that the units module is defined and has the correct version
    const auto [unitsGroup, unitsGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/h5md/modules/units"));
    EXPECT_EQ(getAttributeVector<int>(unitsGroup, "version").value(), (std::vector<int>{ 1, 0 }))
            << "units module specification version should be 1.0";
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

using H5mdSetupFromExistingFile = H5mdIoTest;

TEST_F(H5mdSetupFromExistingFile, WorksForTrajectoryData)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/particles/system"));
    H5mdTimeDataBlockBuilder<RVec>(group, "position").withFrameDimension({ 1 }).build();
    H5mdTimeDataBlockBuilder<RVec>(group, "velocity").withFrameDimension({ 1 }).build();
    H5mdTimeDataBlockBuilder<RVec>(group, "force").withFrameDimension({ 1 }).build();
    const auto [boxGroup, boxGroupGuard] = makeH5mdGroupGuard(createGroup(group, "box/edges"));
    H5mdFrameDataSetBuilder<real>(boxGroup, "value").withFrameDimension({ DIM, DIM }).build();
    EXPECT_NO_THROW(file().setupFromExistingFile());
}

TEST_F(H5mdSetupFromExistingFile, WorksForPositionPlusBoxDataOnly)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/particles/system"));
    H5mdTimeDataBlockBuilder<RVec>(group, "position").withFrameDimension({ 1 }).build();
    const auto [boxGroup, boxGroupGuard] = makeH5mdGroupGuard(createGroup(group, "box/edges"));
    H5mdFrameDataSetBuilder<real>(boxGroup, "value").withFrameDimension({ DIM, DIM }).build();
    EXPECT_NO_THROW(file().setupFromExistingFile());
}

TEST_F(H5mdSetupFromExistingFile, ThrowsForPositionWithoutBoxData)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/particles/system"));
    H5mdTimeDataBlockBuilder<RVec>(group, "position").withFrameDimension({ 1 }).build();
    EXPECT_THROW(file().setupFromExistingFile(), gmx::FileIOError)
            << "Must throw if there is a position but not a box data block";
}

TEST_F(H5mdSetupFromExistingFile, WorksForVelocityDataOnly)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/particles/system"));
    H5mdTimeDataBlockBuilder<RVec>(group, "velocity").withFrameDimension({ 1 }).build();
    EXPECT_NO_THROW(file().setupFromExistingFile());
}

TEST_F(H5mdSetupFromExistingFile, WorksForForceDataOnly)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/particles/system"));
    H5mdTimeDataBlockBuilder<RVec>(group, "force").withFrameDimension({ 1 }).build();
    EXPECT_NO_THROW(file().setupFromExistingFile());
}

TEST_F(H5mdSetupFromExistingFile, ThrowsIfTrajectoryGroupDoesNotExist)
{
    EXPECT_THROW(file().setupFromExistingFile(), gmx::FileIOError)
            << "Must throw before setting up /particles/system";
    makeH5mdGroupGuard(createGroup(fileid(), "/particles/system"));
    EXPECT_NO_THROW(file().setupFromExistingFile());
}

TEST_F(H5mdSetupFromExistingFile, ThrowsIfTrajectoryDataBlocksHaveInconsistentNumParticles)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/particles/system"));
    H5mdTimeDataBlockBuilder<RVec>(group, "position").withFrameDimension({ 1 }).build();
    H5mdTimeDataBlockBuilder<RVec>(group, "velocity").withFrameDimension({ 2 }).build();
    H5mdTimeDataBlockBuilder<RVec>(group, "force").withFrameDimension({ 1 }).build();
    const auto [boxGroup, boxGroupGuard] = makeH5mdGroupGuard(createGroup(group, "box/edges"));
    H5mdFrameDataSetBuilder<real>(boxGroup, "value").withFrameDimension({ DIM, DIM }).build();
    EXPECT_THROW(file().setupFromExistingFile(), gmx::FileIOError)
            << "Must throw if blocks have different numParticles";
}

TEST_F(H5mdIoTest, SetupFileFromInputSetsUnitsToTrajectoryDataSets)
{
    gmx_mtop_t mtop;
    mtop.natoms = 1;
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 1;
    file().setupFileFromInput(mtop, inputRecord);

    const auto [positionGroup, positionGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/position"));
    const auto [velocityGroup, velocityGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/velocity"));
    const auto [forceGroup, forceGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/force"));
    const auto [boxGroup, boxGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid(), "/particles/system/box/edges"));

    // Open the value data sets as raw handles to check their unit attributes
    const auto [position, positionGuard] =
            makeH5mdDataSetGuard(H5Dopen(positionGroup, "value", H5P_DEFAULT));
    EXPECT_EQ(getAttribute<std::string>(position, "unit").value_or(""), "nm");

    const auto [velocity, velocityGuard] =
            makeH5mdDataSetGuard(H5Dopen(velocityGroup, "value", H5P_DEFAULT));
    EXPECT_EQ(getAttribute<std::string>(velocity, "unit").value_or(""), "nm ps-1");

    const auto [force, forceGuard] = makeH5mdDataSetGuard(H5Dopen(forceGroup, "value", H5P_DEFAULT));
    EXPECT_EQ(getAttribute<std::string>(force, "unit").value_or(""), "kJ mol-1 nm-1");

    const auto [box, boxGuard] = makeH5mdDataSetGuard(H5Dopen(boxGroup, "value", H5P_DEFAULT));
    EXPECT_EQ(getAttribute<std::string>(box, "unit").value_or(""), "nm");
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

TEST_F(H5mdIoTest, WriteNextFrameWorks)
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

    constexpr int                            numFrames = 3;
    const std::array<int64_t, numFrames>     steps     = { 5001, 10, 200050 };
    const std::array<double, numFrames>      times     = { 50.1, -1.15, 256.1 };
    std::array<std::vector<RVec>, numFrames> positions;
    std::array<std::vector<RVec>, numFrames> velocities;
    std::array<std::vector<RVec>, numFrames> forces;
    std::array<matrix, numFrames>            boxes;

    // Write unique per-frame values to the file (and to our per-frame buffers above for verification)
    for (int frameIndex = 0; frameIndex < numFrames; ++frameIndex)
    {
        positions[frameIndex].resize(numAtoms);
        velocities[frameIndex].resize(numAtoms);
        forces[frameIndex].resize(numAtoms);
        // For each atom in this frame, create a unique value and set it on all per-frame vectors
        RVec frameValue = static_cast<real>(frameIndex) * RVec{ 1.0, 0.1, 0.01 };
        for (hsize_t atomIndex = 0; atomIndex < numAtoms; ++atomIndex)
        {
            frameValue += { 1.0, 1.0, 1.0 };
            positions[frameIndex][atomIndex]  = frameValue;
            velocities[frameIndex][atomIndex] = static_cast<real>(10.0) * frameValue;
            forces[frameIndex][atomIndex]     = static_cast<real>(100.0) * frameValue;
        }

        // Generate unique values for the box matrix
        for (int i = 0; i < DIM; ++i)
        {
            for (int j = 0; j < DIM; ++j)
            {
                boxes[frameIndex][i][j] = (9 * frameIndex) + (3 * i) + j;
            }
        }

        file().writeNextFrame(positions[frameIndex],
                              velocities[frameIndex],
                              forces[frameIndex],
                              boxes[frameIndex],
                              steps[frameIndex],
                              times[frameIndex]);
    }

    H5mdTimeDataBlock<RVec> positionDataSet(fileid(), "/particles/system/position");
    H5mdTimeDataBlock<RVec> velocityDataSet(fileid(), "/particles/system/velocity");
    H5mdTimeDataBlock<RVec> forceDataSet(fileid(), "/particles/system/force");
    H5mdFrameDataSet<real>  boxDataSet(fileid(), "/particles/system/box/edges/value");

    {
        SCOPED_TRACE("Assert that the correct number of frames were written");
        EXPECT_EQ(positionDataSet.numFrames(), numFrames);
        EXPECT_EQ(velocityDataSet.numFrames(), numFrames);
        EXPECT_EQ(forceDataSet.numFrames(), numFrames);
        EXPECT_EQ(boxDataSet.numFrames(), numFrames);
    }
    {
        SCOPED_TRACE("Assert that trajectory data was written correctly");
        for (int i = 0; i < numFrames; ++i)
        {
            std::vector<RVec> readValueBuffer(numAtoms);
            positionDataSet.readValueAtIndex(i, readValueBuffer);
            EXPECT_EQ(readValueBuffer, positions[i]);
            velocityDataSet.readValueAtIndex(i, readValueBuffer);
            EXPECT_EQ(readValueBuffer, velocities[i]);
            forceDataSet.readValueAtIndex(i, readValueBuffer);
            EXPECT_EQ(readValueBuffer, forces[i]);
            std::array<real, DIM * DIM> readBoxBuffer;
            boxDataSet.readFrame(i, readBoxBuffer);
            EXPECT_THAT(constArrayRefFromArray(reinterpret_cast<real*>(boxes[i]), 9),
                        ::testing::Pointwise(::testing::Eq(), readBoxBuffer));
        }
    }
    {
        SCOPED_TRACE("Assert that steps were written correctly");
        for (int i = 0; i < numFrames; ++i)
        {
            EXPECT_EQ(*positionDataSet.readStepAtIndex(i), steps[i]);
            EXPECT_EQ(*velocityDataSet.readStepAtIndex(i), steps[i]);
            EXPECT_EQ(*forceDataSet.readStepAtIndex(i), steps[i]);
        }
    }
    {
        SCOPED_TRACE("Assert that times were written correctly");
        for (int i = 0; i < numFrames; ++i)
        {
            EXPECT_FLOAT_EQ(*positionDataSet.readTimeAtIndex(i), times[i]);
            EXPECT_FLOAT_EQ(*velocityDataSet.readTimeAtIndex(i), times[i]);
            EXPECT_FLOAT_EQ(*forceDataSet.readTimeAtIndex(i), times[i]);
        }
    }
}

TEST_F(H5mdIoTest, WriteNextFrameDoesNotWriteEmptyRefs)
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

    std::vector<RVec> positionsToWrite(numAtoms);
    std::vector<RVec> velocitiesToWrite(numAtoms);
    std::vector<RVec> forcesToWrite(numAtoms);

    // Generate unique values to fill each value array above
    RVec atomValue = RVec{ 1.0, 0.1, 0.01 };
    for (hsize_t atomIndex = 0; atomIndex < numAtoms; ++atomIndex)
    {
        atomValue += { 1.0, 1.0, 1.0 };
        positionsToWrite[atomIndex]  = atomValue;
        velocitiesToWrite[atomIndex] = static_cast<real>(10.0) * atomValue;
        forcesToWrite[atomIndex]     = static_cast<real>(100.0) * atomValue;
    }

    file().writeNextFrame({}, {}, {}, c_unusedBox, 0, 0);
    file().writeNextFrame(positionsToWrite, {}, {}, c_unusedBox, 0, 0);
    file().writeNextFrame({}, velocitiesToWrite, {}, c_unusedBox, 0, 0);
    file().writeNextFrame({}, {}, forcesToWrite, c_unusedBox, 0, 0);

    H5mdTimeDataBlock<RVec> positionDataSet(fileid(), "/particles/system/position");
    H5mdTimeDataBlock<RVec> velocityDataSet(fileid(), "/particles/system/velocity");
    H5mdTimeDataBlock<RVec> forceDataSet(fileid(), "/particles/system/force");
    {
        SCOPED_TRACE("Assert that only one frame was written to each data set");
        EXPECT_EQ(positionDataSet.numFrames(), 1);
        EXPECT_EQ(velocityDataSet.numFrames(), 1);
        EXPECT_EQ(forceDataSet.numFrames(), 1);
    }
    {
        SCOPED_TRACE("Assert that the correct values were written to frame 0");
        std::vector<RVec> readValueBuffer(numAtoms);
        positionDataSet.readValueAtIndex(0, readValueBuffer);
        EXPECT_EQ(readValueBuffer, positionsToWrite);
        velocityDataSet.readValueAtIndex(0, readValueBuffer);
        EXPECT_EQ(readValueBuffer, velocitiesToWrite);
        forceDataSet.readValueAtIndex(0, readValueBuffer);
        EXPECT_EQ(readValueBuffer, forcesToWrite);
    }
}

TEST_F(H5mdIoTest, WriteNextFrameThrowsForNotCreatedDataSets)
{
    t_inputrec inputRecord;
    // Do not construct data sets for any trajectory data
    inputRecord.nstxout = 0;
    inputRecord.nstvout = 0;
    inputRecord.nstfout = 0;

    // Read the topology from a test system in our simulation data base
    const std::string fileNameBase = "spc2-traj";
    TprAndFileManager tprFileHandle(fileNameBase);
    bool              haveTopology;
    gmx_mtop_t        mtop;
    readConfAndTopology(tprFileHandle.tprName(), &haveTopology, &mtop, nullptr, nullptr, nullptr, nullptr);
    const hsize_t numAtoms = static_cast<hsize_t>(mtop.natoms);

    file().setupFileFromInput(mtop, inputRecord);

    std::vector<RVec> valuesToWrite(numAtoms, { 0.0, 0.0, 0.0 });
    ASSERT_NO_THROW(file().writeNextFrame({}, {}, {}, c_unusedBox, 0, 0))
            << "Sanity check failed: should not throw when not trying to write any data";
    EXPECT_THROW(file().writeNextFrame(valuesToWrite, {}, {}, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_THROW(file().writeNextFrame({}, valuesToWrite, {}, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_THROW(file().writeNextFrame({}, {}, valuesToWrite, c_unusedBox, 0, 0), gmx::FileIOError);
}

TEST_F(H5mdIoTest, WriteNextFrameThrowsForBuffersWithIncorrectSize)
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

    std::vector<RVec> bufferTooSmall(numAtoms - 1, { 0.0, 0.0, 0.0 });
    std::vector<RVec> bufferTooLarge(numAtoms + 1, { 0.0, 0.0, 0.0 });
    std::vector<RVec> bufferJustRight(numAtoms, { 0.0, 0.0, 0.0 });

    EXPECT_THROW(file().writeNextFrame(bufferTooSmall, {}, {}, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_THROW(file().writeNextFrame(bufferTooLarge, {}, {}, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_NO_THROW(file().writeNextFrame(bufferJustRight, {}, {}, c_unusedBox, 0, 0));
    EXPECT_THROW(file().writeNextFrame({}, bufferTooSmall, {}, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_THROW(file().writeNextFrame({}, bufferTooLarge, {}, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_NO_THROW(file().writeNextFrame({}, bufferJustRight, {}, c_unusedBox, 0, 0));
    EXPECT_THROW(file().writeNextFrame({}, {}, bufferTooSmall, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_THROW(file().writeNextFrame({}, {}, bufferTooLarge, c_unusedBox, 0, 0), gmx::FileIOError);
    EXPECT_NO_THROW(file().writeNextFrame({}, {}, bufferJustRight, c_unusedBox, 0, 0));
}

//! \brief Helper function to return an ArrayRef<RVec> from an input rvec pointer \p values.
ArrayRef<RVec> asRVecArray(rvec* values, int64_t numValues)
{
    return arrayRefFromArray(reinterpret_cast<RVec*>(values), numValues);
}

using H5mdReadNextFrame = H5mdIoTest;

TEST_F(H5mdReadNextFrame, Works)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 1;

    const int  numAtoms = 6;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);

    constexpr int                            numFrames = 3;
    const std::array<int64_t, numFrames>     steps     = { 5001, 10, 200050 };
    const std::array<double, numFrames>      times     = { 50.1, -1.15, 256.1 };
    std::array<std::vector<RVec>, numFrames> positions;
    std::array<std::vector<RVec>, numFrames> velocities;
    std::array<std::vector<RVec>, numFrames> forces;
    std::array<matrix, numFrames>            boxes;

    // Write unique per-frame values to the file (and to our per-frame buffers above for verification)
    for (int frameIndex = 0; frameIndex < numFrames; ++frameIndex)
    {
        positions[frameIndex].resize(numAtoms);
        velocities[frameIndex].resize(numAtoms);
        forces[frameIndex].resize(numAtoms);
        // For each atom in this frame, create a unique value and append to all per-frame vectors
        RVec frameValue = static_cast<real>(frameIndex) * RVec{ 1.0, 0.1, 0.01 };
        for (hsize_t atomIndex = 0; atomIndex < numAtoms; ++atomIndex)
        {
            frameValue += { 1.0, 1.0, 1.0 };
            positions[frameIndex][atomIndex]  = frameValue;
            velocities[frameIndex][atomIndex] = static_cast<real>(10.0) * frameValue;
            forces[frameIndex][atomIndex]     = static_cast<real>(100.0) * frameValue;
        }

        // Generate unique values for the box matrix
        for (int i = 0; i < DIM; ++i)
        {
            for (int j = 0; j < DIM; ++j)
            {
                boxes[frameIndex][i][j] = (9 * frameIndex) + (3 * i) + j;
            }
        }


        file().writeNextFrame(positions[frameIndex],
                              velocities[frameIndex],
                              forces[frameIndex],
                              boxes[frameIndex],
                              steps[frameIndex],
                              times[frameIndex]);
    }

    t_trxframe* frame;
    snew(frame, 1);
    for (int frameIndex = 0; frameIndex < numFrames; ++frameIndex)
    {
        EXPECT_TRUE(file().readNextFrame(frame)) << "Must return true for a successful read";

        EXPECT_TRUE(frame->bX);
        EXPECT_THAT(asRVecArray(frame->x, frame->natoms),
                    ::testing::Pointwise(::testing::Eq(), positions[frameIndex]));
        EXPECT_TRUE(frame->bV);
        EXPECT_THAT(asRVecArray(frame->v, frame->natoms),
                    ::testing::Pointwise(::testing::Eq(), velocities[frameIndex]));
        EXPECT_TRUE(frame->bF);
        EXPECT_THAT(asRVecArray(frame->f, frame->natoms),
                    ::testing::Pointwise(::testing::Eq(), forces[frameIndex]));
        EXPECT_TRUE(frame->bX);
        EXPECT_THAT(constArrayRefFromArray(reinterpret_cast<real*>(frame->box), DIM * DIM),
                    ::testing::Pointwise(
                            ::testing::FloatEq(),
                            constArrayRefFromArray(reinterpret_cast<real*>(boxes[frameIndex]), DIM * DIM)));
        EXPECT_TRUE(frame->bStep);
        EXPECT_EQ(frame->step, steps[frameIndex]);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, times[frameIndex]);
    }
    EXPECT_FALSE(file().readNextFrame(frame))
            << "Must return false when no more frames exist to read";
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, ReturnsFalseBeforeDataIsWritten)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 1;
    const int  numAtoms = 1;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;
    file().setupFileFromInput(mtop, inputRecord);

    t_trxframe* frame;
    snew(frame, 1);
    EXPECT_FALSE(file().readNextFrame(frame));
    std::vector<RVec> valuesToWrite = { { 0.0, 1.0, 2.0 } };
    file().writeNextFrame(valuesToWrite, valuesToWrite, valuesToWrite, c_unusedBox, 0, 0.0);
    EXPECT_TRUE(file().readNextFrame(frame));
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, WorksIfNoDataSetsExists)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 0;
    inputRecord.nstvout = 0;
    inputRecord.nstfout = 0;
    const int  numAtoms = 1;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;
    file().setupFileFromInput(mtop, inputRecord);

    t_trxframe* frame;
    snew(frame, 1);
    EXPECT_FALSE(file().readNextFrame(frame));
    EXPECT_FALSE(frame->bX);
    EXPECT_FALSE(frame->bV);
    EXPECT_FALSE(frame->bF);
    EXPECT_FALSE(frame->bBox);
    EXPECT_FALSE(frame->bStep);
    EXPECT_FALSE(frame->bTime);
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, WorksIfOnlyPositionDataExists)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 0;
    inputRecord.nstfout = 0;
    const int  numAtoms = 1;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);

    const std::vector<RVec> valuesToWrite = { { 0.0, 1.0, 2.0 } };
    matrix boxToWrite = { { 0.1, 0.2, 0.3 }, { 1.1, 1.2, 1.3 }, { 2.1, 2.2, 2.3 } };
    file().writeNextFrame(valuesToWrite, {}, {}, boxToWrite, -1, -1.0);

    t_trxframe* frame;
    snew(frame, 1);
    EXPECT_TRUE(file().readNextFrame(frame));
    EXPECT_TRUE(frame->bX);
    EXPECT_THAT(asRVecArray(frame->x, frame->natoms), ::testing::Pointwise(::testing::Eq(), valuesToWrite));
    EXPECT_FALSE(frame->bV);
    EXPECT_FALSE(frame->bF);
    EXPECT_TRUE(frame->bBox) << "Simulation box must be read when positions are";
    EXPECT_THAT(constArrayRefFromArray(reinterpret_cast<real*>(frame->box), DIM * DIM),
                ::testing::Pointwise(
                        ::testing::FloatEq(),
                        constArrayRefFromArray(reinterpret_cast<real*>(boxToWrite), DIM * DIM)));
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, WorksIfOnlyVelocityDataExists)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 0;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 0;
    const int  numAtoms = 1;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);

    const std::vector<RVec> valuesToWrite = { { 0.0, 1.0, 2.0 } };
    file().writeNextFrame({}, valuesToWrite, {}, c_unusedBox, -1, -1.0);

    t_trxframe* frame;
    snew(frame, 1);
    EXPECT_TRUE(file().readNextFrame(frame));
    EXPECT_FALSE(frame->bX);
    EXPECT_TRUE(frame->bV);
    EXPECT_THAT(asRVecArray(frame->v, frame->natoms), ::testing::Pointwise(::testing::Eq(), valuesToWrite));
    EXPECT_FALSE(frame->bF);
    EXPECT_FALSE(frame->bBox);
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, WorksIfOnlyForceDataExists)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 0;
    inputRecord.nstvout = 0;
    inputRecord.nstfout = 1;
    const int  numAtoms = 1;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;

    file().setupFileFromInput(mtop, inputRecord);

    const std::vector<RVec> valuesToWrite = { { 0.0, 1.0, 2.0 } };
    matrix boxToWrite = { { 0.1, 0.2, 0.3 }, { 1.1, 1.2, 1.3 }, { 2.1, 2.2, 2.3 } };
    file().writeNextFrame({}, {}, valuesToWrite, boxToWrite, -1, -1.0);

    t_trxframe* frame;
    snew(frame, 1);
    EXPECT_TRUE(file().readNextFrame(frame));
    EXPECT_FALSE(frame->bX);
    EXPECT_FALSE(frame->bV);
    EXPECT_TRUE(frame->bF);
    EXPECT_THAT(asRVecArray(frame->f, frame->natoms), ::testing::Pointwise(::testing::Eq(), valuesToWrite));
    EXPECT_FALSE(frame->bBox);
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, DataSetsWithDifferentStepFrequenciesAreReadInOrder)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 1;

    const int  numAtoms = 1;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;
    file().setupFileFromInput(mtop, inputRecord);

    // First frame: Velocity only
    const int64_t           step1     = 1;
    const std::vector<RVec> velocity1 = { { 0.3, 0.4, 0.5 } };
    file().writeNextFrame({}, velocity1, {}, c_unusedBox, step1, static_cast<double>(step1));

    // Second frame: Position + Force
    const int64_t           step2     = 2;
    const std::vector<RVec> position1 = { { 0.0, 0.1, 0.2 } };
    const std::vector<RVec> force1    = { { 0.6, 0.7, 0.8 } };
    file().writeNextFrame(position1, {}, force1, c_unusedBox, step2, static_cast<double>(step2));

    // Third frame: Velocity + force
    const int64_t           step3     = 3;
    const std::vector<RVec> velocity2 = { { 1.3, 1.4, 1.5 } };
    const std::vector<RVec> force2    = { { 1.6, 1.7, 1.8 } };
    file().writeNextFrame({}, velocity2, force2, c_unusedBox, step3, static_cast<double>(step3));

    // Fourth frame: Position only
    const int64_t           step4     = 4;
    const std::vector<RVec> position2 = { { 1.0, 1.1, 1.2 } };
    file().writeNextFrame(position2, {}, {}, c_unusedBox, step4, static_cast<double>(step4));

    // Fifth frame: All
    const int64_t           step5     = 5;
    const std::vector<RVec> position3 = { { 2.0, 2.1, 2.2 } };
    const std::vector<RVec> velocity3 = { { 2.3, 2.4, 2.5 } };
    const std::vector<RVec> force3    = { { 2.6, 2.7, 2.8 } };
    file().writeNextFrame(position3, velocity3, force3, c_unusedBox, step5, static_cast<double>(step5));

    t_trxframe* frame;
    snew(frame, 1);
    {
        SCOPED_TRACE("Frame 1: Velocity only");
        file().readNextFrame(frame);
        EXPECT_FALSE(frame->bX);
        EXPECT_TRUE(frame->bV);
        EXPECT_THAT(asRVecArray(frame->v, frame->natoms), ::testing::Pointwise(::testing::Eq(), velocity1));
        EXPECT_FALSE(frame->bF);
        EXPECT_FALSE(frame->bBox) << "Box frames are only read when positions are read";
        EXPECT_TRUE(frame->bStep);
        EXPECT_EQ(frame->step, step1);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, static_cast<double>(step1));
    }
    {
        SCOPED_TRACE("Frame 2: Position + Force");
        file().readNextFrame(frame);
        EXPECT_TRUE(frame->bX);
        EXPECT_THAT(asRVecArray(frame->x, frame->natoms), ::testing::Pointwise(::testing::Eq(), position1));
        EXPECT_FALSE(frame->bV);
        EXPECT_TRUE(frame->bF);
        EXPECT_THAT(asRVecArray(frame->f, frame->natoms), ::testing::Pointwise(::testing::Eq(), force1));
        EXPECT_TRUE(frame->bBox) << "Box frames are only read when positions are read";
        EXPECT_TRUE(frame->bStep);
        EXPECT_EQ(frame->step, step2);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, static_cast<double>(step2));
    }
    {
        SCOPED_TRACE("Frame 3: Velocity + Force");
        file().readNextFrame(frame);
        EXPECT_FALSE(frame->bX);
        EXPECT_TRUE(frame->bV);
        EXPECT_THAT(asRVecArray(frame->v, frame->natoms), ::testing::Pointwise(::testing::Eq(), velocity2));
        EXPECT_TRUE(frame->bF);
        EXPECT_THAT(asRVecArray(frame->f, frame->natoms), ::testing::Pointwise(::testing::Eq(), force2));
        EXPECT_FALSE(frame->bBox) << "Box frames are only read when positions are read";
        EXPECT_TRUE(frame->bStep);
        EXPECT_EQ(frame->step, step3);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, static_cast<double>(step3));
    }
    {
        SCOPED_TRACE("Frame 4: Position only");
        file().readNextFrame(frame);
        EXPECT_TRUE(frame->bX);
        EXPECT_THAT(asRVecArray(frame->x, frame->natoms), ::testing::Pointwise(::testing::Eq(), position2));
        EXPECT_FALSE(frame->bV);
        EXPECT_FALSE(frame->bF);
        EXPECT_TRUE(frame->bBox) << "Box frames are only read when positions are read";
        EXPECT_TRUE(frame->bStep);
        EXPECT_EQ(frame->step, step4);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, static_cast<double>(step4));
    }
    {
        SCOPED_TRACE("Frame 3: Position + Velocity + Force");
        file().readNextFrame(frame);
        EXPECT_TRUE(frame->bX);
        EXPECT_THAT(asRVecArray(frame->x, frame->natoms), ::testing::Pointwise(::testing::Eq(), position3));
        EXPECT_TRUE(frame->bV);
        EXPECT_THAT(asRVecArray(frame->v, frame->natoms), ::testing::Pointwise(::testing::Eq(), velocity3));
        EXPECT_TRUE(frame->bF);
        EXPECT_THAT(asRVecArray(frame->f, frame->natoms), ::testing::Pointwise(::testing::Eq(), force3));
        EXPECT_TRUE(frame->bBox) << "Box frames are only read when positions are read";
        EXPECT_TRUE(frame->bStep);
        EXPECT_EQ(frame->step, step5);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, static_cast<double>(step5));
    }
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, MissingTimeDataSetsAreHandled)
{
    // Manually set up the trajectory tree, only position and force data sets store time!
    const DataSetDims frameDims = { 1 };
    {
        const auto [positionGroup, positionGroupGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/particles/system/position"));
        H5mdFrameDataSetBuilder<RVec>(positionGroup, "value").withFrameDimension(frameDims).build();
        H5mdFrameDataSetBuilder<int64_t>(positionGroup, "step").build();
        H5mdFrameDataSetBuilder<double>(positionGroup, "time").build();

        const auto [velocityGroup, velocityGroupGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/particles/system/velocity"));
        H5mdFrameDataSetBuilder<RVec>(velocityGroup, "value").withFrameDimension(frameDims).build();
        H5mdFrameDataSetBuilder<int64_t>(velocityGroup, "step").build();
        // No time data set here!

        const auto [forceGroup, forceGroupGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/particles/system/force"));
        H5mdFrameDataSetBuilder<RVec>(forceGroup, "value").withFrameDimension(frameDims).build();
        H5mdFrameDataSetBuilder<int64_t>(forceGroup, "step").build();
        H5mdFrameDataSetBuilder<double>(forceGroup, "time").build();

        const auto [boxGroup, boxGroupGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/particles/system/box/edges"));
        H5mdFrameDataSetBuilder<real>(boxGroup, "value").withFrameDimension({ DIM, DIM }).build();

        file().setupFromExistingFile();
    }

    const std::vector<RVec> valuesToWrite = { { -1.0, -1.0, -1.0 } };

    // Frame 1: Write velocities only (will not write time)
    const double time1 = 1.1;
    file().writeNextFrame({}, valuesToWrite, {}, c_unusedBox, 1, time1);

    // Frame 2: Write position + velocity (will write time)
    const double time2 = 2.2;
    file().writeNextFrame(valuesToWrite, valuesToWrite, {}, c_unusedBox, 2, time2);

    // Frame 3: Write position only (will write time)
    const double time3 = 3.3;
    file().writeNextFrame(valuesToWrite, {}, {}, c_unusedBox, 3, time3);

    // Frame 4: Write velocities only (will not write time)
    const double time4 = 4.4;
    file().writeNextFrame({}, valuesToWrite, {}, c_unusedBox, 4, time4);

    // Frame 5: Write velocities + force (will write time)
    const double time5 = 5.5;
    file().writeNextFrame({}, valuesToWrite, valuesToWrite, c_unusedBox, 5, time5);

    // Frame 1: Write force only (will write time)
    const double time6 = 6.6;
    file().writeNextFrame({}, {}, valuesToWrite, c_unusedBox, 6, time6);

    t_trxframe* frame;
    snew(frame, 1);
    {
        SCOPED_TRACE("Read frame 1 (velocity only, will not have time)");
        file().readNextFrame(frame);
        EXPECT_FALSE(frame->bTime);
    }
    {
        SCOPED_TRACE("Read frame 2 (position + velocity, will have time)");
        file().readNextFrame(frame);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, time2);
    }
    {
        SCOPED_TRACE("Read frame 3 (position only, will have time)");
        file().readNextFrame(frame);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, time3);
    }
    {
        SCOPED_TRACE("Read frame 4 (velocity only, will not have time)");
        file().readNextFrame(frame);
        EXPECT_FALSE(frame->bTime);
    }
    {
        SCOPED_TRACE("Read frame 5 (velocity + force, will have time)");
        file().readNextFrame(frame);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, time5);
    }
    {
        SCOPED_TRACE("Read frame 6 (force only, will have time)");
        file().readNextFrame(frame);
        EXPECT_TRUE(frame->bTime);
        EXPECT_FLOAT_EQ(frame->time, time6);
    }
    done_frame(frame);
}

TEST_F(H5mdReadNextFrame, NonTrajectoryFrameBoolsInTrxFrameAreFalse)
{
    t_inputrec inputRecord;
    inputRecord.nstxout = 1;
    inputRecord.nstvout = 1;
    inputRecord.nstfout = 1;

    const int  numAtoms = 1;
    gmx_mtop_t mtop;
    mtop.natoms = numAtoms;
    file().setupFileFromInput(mtop, inputRecord);

    const std::vector<RVec> valuesToWrite = { { -1.0, -1.0, -1.0 } };
    file().writeNextFrame(valuesToWrite, valuesToWrite, valuesToWrite, c_unusedBox, -1, -1.0);

    t_trxframe* frame;
    snew(frame, 1);

    file().readNextFrame(frame);
    EXPECT_FALSE(frame->bAtoms);
    EXPECT_FALSE(frame->bLambda);
    EXPECT_FALSE(frame->bFepState);
    EXPECT_FALSE(frame->bIndex);
    EXPECT_FALSE(frame->bPBC);
    EXPECT_FALSE(frame->bPrec); // check this for reduced-precision trajectories separately
#if GMX_DOUBLE
    EXPECT_TRUE(frame->bDouble);
#else
    EXPECT_FALSE(frame->bDouble);
#endif

    done_frame(frame);
}

} // namespace
} // namespace test
} // namespace gmx
