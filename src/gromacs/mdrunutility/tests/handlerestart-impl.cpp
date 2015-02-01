/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements unit tests for handling mdrun restarts
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "gromacs/mdrunutility/handlerestart-impl.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

#include "mockfileiohandler.h"
#include "pargswrapper.h"

namespace gmx
{

namespace HandleRestart
{

namespace test
{

using gmx::test::MockFileIOHandler;
using testing::_;
using testing::Return;
using testing::SetArgPointee;
using testing::StrEq;

//! Base class for test fixtures that test some methods of Impl
class ImplTestBase : public ::testing::Test
{
    public:
        ImplTestBase() : input_(),
                         defaultData_()
        {
        }

        //! Mock input object
        MockFileIOHandler input_;
        /*! \brief Default-constructed object of data read from
            checkpoint file, used in test assertions */
        DataFromCheckpoint defaultData_;
};

/* Ideally, these would be static const members of
   GetDataFromCheckpointTest, but C++98 can't cope with in-class array
   initializers, which means that the ArrayRef template constructors
   that are based on these arrays can't work. */
//! Fake mdrun argv array, used in most tests
const char *GetDataFromCheckpointTestArgv_g[] = { "mdrun", "-cpi", "state.cpt", "-g", "md.log" };
//! Fake mdrun filenm structure to fill from argv, used in most tests
t_filenm    GetDataFromCheckpointTestFilenames_g[] = { { efCPT, "-cpi", NULL, ffSET },
                                                       { efLOG, "-g", NULL, ffSET } };

//! Test fixture for getDataFromCheckpointTest() method of Impl
class GetDataFromCheckpointTest : public ImplTestBase
{
    public:
        GetDataFromCheckpointTest() : ImplTestBase(),
                                      impl_(&input_, NULL, ConstCommandLineFilenames(GetDataFromCheckpointTestFilenames_g))
        {
        }
        //! Called once for all tests in test case
        static void SetUpTestCase()
        {
            /* Build the rest of the filenames array from the argv
               inputs. Only needs to be done once for all test cases,
               because the data is logically const. */
            parseCommonArgs(GetDataFromCheckpointTestArgv_g, GetDataFromCheckpointTestFilenames_g);
        }
        //! Impl object normally used in the tests
        Impl impl_;
};

TEST_F(GetDataFromCheckpointTest, ReturnsDefaultDataIfNoCpiOptionWasGiven)
{
    // Get the set of filenames, but skip the -cpi option, which is first
    ConstCommandLineFilenames filenamesWithoutCpiOption = ConstCommandLineFilenames(GetDataFromCheckpointTestFilenames_g+1,
                                                                                    GetDataFromCheckpointTestFilenames_g+2);
    Impl impl(&input_, NULL, filenamesWithoutCpiOption);

    // Note, not using impl_
    DataFromCheckpoint actualData = impl.getDataFromCheckpoint(true, false);
    // Check that default values were returned
    EXPECT_EQ(defaultData_.partNumber_, actualData.partNumber_);
}

TEST_F(GetDataFromCheckpointTest, ReturnsDefaultIfNotSimMaster)
{
    // Use the normal filenames in impl_
    bool               bIsSimMaster = false;
    DataFromCheckpoint actualData   = impl_.getDataFromCheckpoint(bIsSimMaster, false);
    // Check that default values were returned
    EXPECT_EQ(defaultData_.partNumber_, actualData.partNumber_);
}

TEST_F(GetDataFromCheckpointTest, ReturnsDefaultIfNeitherCheckpointOrLogFileExist)
{
    // Use the normal filenames in impl_
    EXPECT_CALL(input_, gmx_fexist(StrEq("state.cpt"))).WillOnce(Return(false));
    EXPECT_CALL(input_, gmx_fexist(StrEq("md.log"))).WillOnce(Return(false));
    DataFromCheckpoint actualData = impl_.getDataFromCheckpoint(true, false);
    // Check that default values were returned
    EXPECT_EQ(defaultData_.partNumber_, actualData.partNumber_);
}

TEST_F(GetDataFromCheckpointTest, ThrowsIfCheckpointFileDoesNotExistButLogFileDoes)
{
    // Use the normal filenames in impl_
    EXPECT_CALL(input_, gmx_fexist(StrEq("state.cpt"))).WillOnce(Return(false));
    EXPECT_CALL(input_, gmx_fexist(StrEq("md.log"))).WillOnce(Return(true));
    EXPECT_THROW_GMX(impl_.getDataFromCheckpoint(true, false), InconsistentInputError);
}

TEST_F(GetDataFromCheckpointTest, ThrowsWhenCheckpointHadNoOutputFilenames)
{
    // Use the normal filenames in impl_
    int                  numFilenamesInCheckpoint = 0;
    gmx_file_position_t *empty                    = NULL;

    EXPECT_CALL(input_, gmx_fexist(StrEq("state.cpt"))).WillOnce(Return(true));
    EXPECT_CALL(input_, read_checkpoint_simulation_part_and_filenames(StrEq("state.cpt"), _, _, _))
        .WillOnce(DoAll(SetArgPointee<2>(numFilenamesInCheckpoint),
                        SetArgPointee<3>(empty)));
    EXPECT_THROW_GMX(impl_.getDataFromCheckpoint(true, false), InternalError);
}

TEST_F(GetDataFromCheckpointTest, ReadsCheckpointDataWhenAppropriate)
{
    // Use the normal filenames in impl_
    int                  checkpointPart = 10;
    gmx_file_position_t *empty          = NULL;

    EXPECT_CALL(input_, gmx_fexist(StrEq("state.cpt"))).WillOnce(Return(true));
    EXPECT_CALL(input_, read_checkpoint_simulation_part_and_filenames(StrEq("state.cpt"), _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(-1),
                        SetArgPointee<3>(empty)));
    DataFromCheckpoint actualData = impl_.getDataFromCheckpoint(true, false);
    // Check that default values were returned
    EXPECT_EQ(checkpointPart, actualData.partNumber_);
}

//! Fake mdrun argv array, used in most tests
const char *OutputFileExistsTestArgv_g[] = { "mdrun", "-o1", "firstOutput.log", "-o2", "secondOutput.xvg" };
//! Fake mdrun filenm structure to fill from argv, used in most tests
t_filenm    OutputFileExistsTestFilenames_g[] = {
    { efLOG, "-o1",  NULL, ffOPTWR },
    { efXVG, "-o2",  NULL, ffOPTWR },
};

//! Test fixture for outputFileExists() method of Impl
class OutputFileExistsTest : public ImplTestBase
{
    public:
        OutputFileExistsTest() : ImplTestBase(),
                                 impl_(&input_, NULL, ConstCommandLineFilenames(OutputFileExistsTestFilenames_g))
        {
        }
        //! Called once for all tests in test case
        static void SetUpTestCase()
        {
            /* Build the rest of the filenames array from the argv
               inputs. Only needs to be done once for all test cases,
               because the data is logically const. */
            parseCommonArgs(OutputFileExistsTestArgv_g, OutputFileExistsTestFilenames_g);
        }
        //! Impl object normally used in the tests
        Impl impl_;
};

TEST_F(OutputFileExistsTest, CanFindOutputFilesThatExist)
{
    // Use the normal filenames in impl_
    EXPECT_CALL(input_, gmx_fexist(StrEq("firstOutput.log"))).WillOnce(Return(true));
    EXPECT_CALL(input_, gmx_fexist(StrEq("secondOutput.xvg"))).WillOnce(Return(true));
    EXPECT_TRUE(impl_.outputFileExists("firstOutput.log"));
    EXPECT_TRUE(impl_.outputFileExists("secondOutput.xvg"));
}

TEST_F(OutputFileExistsTest, CantFindOutputFilesThatDontExist)
{
    // Use the normal filenames in impl_
    EXPECT_CALL(input_, gmx_fexist(StrEq("secondOutput.xvg"))).WillOnce(Return(false));
    EXPECT_FALSE(impl_.outputFileExists("secondOutput.xvg"));
}

TEST_F(OutputFileExistsTest, CantFindFilesThatAreNotInFnm)
{
    // Use the normal filenames in impl_
    EXPECT_FALSE(impl_.outputFileExists("missing"));
}

/* We can't test that we can't find input files, because a named input
   file in argv & filenames must exist on the filesystem in the normal
   operation of parseCommonArgs. However, this is not an important
   case to test. */
//TEST_F(OutputFileExistsTest, CantFindInputFiles)

TEST_F(OutputFileExistsTest, CantFindFilesInEmptyList)
{
    // Note, not using normal filenames in impl_
    MockFileIOHandler input_;
    Impl              impl(&input_, NULL, EmptyArrayRef());

    EXPECT_FALSE(impl.outputFileExists("inputFilename.xvg"));
}

//! Test fixture for Impl::checkAllFilesForAppendingExist
class CheckAllFilesForAppendingExistTest : public ::testing::Test
{
    public:
        class MockImpl : public Impl
        {
            public:
                MockImpl(FileIOHandlerInterface *input) :
                    Impl(input, NULL, EmptyArrayRef())
                {
                }

                MOCK_CONST_METHOD1(outputFileExists, bool(const std::string&));
                MOCK_CONST_METHOD2(makeMessageWhenSetsOfFilesDontMatch,
                                   std::string(const std::vector<std::string> &,
                                               size_t));

        };
        CheckAllFilesForAppendingExistTest() : impl_(NULL)
        {
            // The checkpoint file must name a log file
            expectedOutputFilenames_.push_back("md.log");
        }

        //! Object under test
        MockImpl                 impl_;
        //! Input data for Impl::checkAllFilesForAppendingExist()
        std::vector<std::string> expectedOutputFilenames_;
};

TEST_F(CheckAllFilesForAppendingExistTest, WorksWhenOutputFilesExist)
{
    EXPECT_CALL(impl_, outputFileExists(StrEq("md.log"))).WillOnce(Return(true));
    EXPECT_NO_THROW(impl_.checkAllFilesForAppendingExist(expectedOutputFilenames_));
}

TEST_F(CheckAllFilesForAppendingExistTest, ThrowsWhenNoOutputFilesExist)
{
    EXPECT_CALL(impl_, outputFileExists(StrEq("md.log"))).WillOnce(Return(false));
    EXPECT_THROW_GMX(impl_.checkAllFilesForAppendingExist(expectedOutputFilenames_),
                     InconsistentInputError);
}

TEST_F(CheckAllFilesForAppendingExistTest, ThrowsWhenNotAllOutputFilesExist)
{
    // Add an output filename that will not exist
    expectedOutputFilenames_.push_back("ener.edr");
    EXPECT_CALL(impl_, outputFileExists(StrEq("md.log"))).WillOnce(Return(true));
    EXPECT_CALL(impl_, outputFileExists(StrEq("ener.edr"))).WillOnce(Return(false));
    EXPECT_CALL(impl_, makeMessageWhenSetsOfFilesDontMatch(expectedOutputFilenames_, 1)).WillOnce(Return(""));
    EXPECT_THROW_GMX(impl_.checkAllFilesForAppendingExist(expectedOutputFilenames_),
                     InconsistentInputError);
}

//! Test fixture for Impl::canAppend
class CanAppendTest : public ::testing::Test
{
    public:
        class MockImpl : public Impl
        {
            public:
                MockImpl(FileIOHandlerInterface *input) :
                    Impl(input, NULL, EmptyArrayRef())
                {
                }

                MOCK_CONST_METHOD1(checkAllFilesForAppendingExist,
                                   void(const std::vector<std::string> &expectedOutputFilenames));
        };
        CanAppendTest() : impl_(NULL)
        {
        }

        //! Object under test
        MockImpl           impl_;
        //! Input data for Impl::canAppend()
        DataFromCheckpoint dataFromCheckpoint_;
};

TEST_F(CanAppendTest, AppendsWhenItShould)
{
    dataFromCheckpoint_.partNumber_ = 1;
    dataFromCheckpoint_.expectedOutputFilenames_.push_back("md.log");
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(dataFromCheckpoint_.expectedOutputFilenames_)).WillOnce(Return());
    EXPECT_TRUE(impl_.canAppend(dataFromCheckpoint_, true));
}

TEST_F(CanAppendTest, WontAppendIfThisIsTheFirstPart)
{
    dataFromCheckpoint_.partNumber_ = 0;
    EXPECT_FALSE(impl_.canAppend(dataFromCheckpoint_, true));
}

TEST_F(CanAppendTest, WontAppendIfShouldNotTry)
{
    dataFromCheckpoint_.partNumber_ = 1;
    bool bTryToAppendFiles = false;
    EXPECT_FALSE(impl_.canAppend(dataFromCheckpoint_, bTryToAppendFiles));
}

TEST_F(CanAppendTest, WontAppendIfCheckpointedSimulationUsedPartNumbers)
{
    dataFromCheckpoint_.partNumber_ = 1;
    dataFromCheckpoint_.expectedOutputFilenames_.push_back("md.part0001.log");
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(dataFromCheckpoint_.expectedOutputFilenames_)).WillOnce(Return());
    EXPECT_FALSE(impl_.canAppend(dataFromCheckpoint_, true));
}

TEST_F(CanAppendTest, ThrowsWhenCheckpointIsMalformed)
{
    dataFromCheckpoint_.partNumber_ = 1;
    /* Set up sets of filenames (ostensibly from the checkpoint and
       command line) that won't match */
    dataFromCheckpoint_.expectedOutputFilenames_.push_back("ener.edr");
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(dataFromCheckpoint_.expectedOutputFilenames_)).WillOnce(Return());
    EXPECT_THROW_GMX(impl_.canAppend(dataFromCheckpoint_, true),
                     InternalError);
}

} // namespace

} // namespace

} // namespace
