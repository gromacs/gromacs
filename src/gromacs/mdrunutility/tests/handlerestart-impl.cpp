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

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

#ifndef DOXYGEN
/* Stub definition, because we need to be able to make
 * one of these during testing */
struct t_fileio {
};
#endif /* DOXYGEN */

namespace gmx
{

namespace HandleRestart
{

namespace test
{

using testing::_;
using testing::Return;
using testing::DoAll;
using testing::SetArgPointee;
using testing::Throw;

/*! \internal
 * \brief Implements the same interface as InputHandler
 *
 * Used when testing Impl::getDataFromCheckpoint() and
 * Impl::outputFileExists()
 */
class TestInputHandler : public InputHandler
{
    public:
        virtual ~TestInputHandler() {}

        MOCK_CONST_METHOD1(gmx_fexist, gmx_bool(const char *));

        MOCK_CONST_METHOD2(gmx_fio_open, t_fileio *(const char *, const char *));

        MOCK_CONST_METHOD4(read_checkpoint_simulation_part_and_filenames,
                           void(t_fileio             *,
                                int                  *,
                                int                  *,
                                gmx_file_position_t **));

        MOCK_CONST_METHOD3(opt2bSet, gmx_bool(const char *, int, const t_filenm []));

        MOCK_CONST_METHOD3(ftp2fn, const char *(int, int, const t_filenm []));

        MOCK_CONST_METHOD3(opt2fn,
                           const char *(const char *, int, const t_filenm []));
};

class TestRankHandler : public RankHandler
{
    public:
        MOCK_CONST_METHOD0(hasMultipleRanks, bool());
        MOCK_CONST_METHOD0(isMaster, bool());
        MOCK_CONST_METHOD0(isSimMaster, bool());
        MOCK_CONST_METHOD0(isMultiSim, bool());
        MOCK_CONST_METHOD0(isMultiMaster, bool());
        MOCK_CONST_METHOD2(broadcast, void(int, void *));
        MOCK_CONST_METHOD4(checkAcrossMultiSim, void(FILE *, int, const char *, bool));
};

//! Test fixture for some methods of Impl
class ImplTest : public ::testing::Test
{
    public:
        ImplTest() : input_(new TestInputHandler),
                     impl_(input_, NULL, EmptyArrayRef()),
                     defaultData_()
        {
        }

        ~ImplTest()
        {
        }

        //! Keep track of the mock input object; impl_ owns it
        TestInputHandler  *input_;
        //! Object under test
        Impl               impl_;
        /*! \brief Default-constructed object of data read from
            checkpoint file, used in test assertions */
        DataFromCheckpoint defaultData_;
};

//! Typedef for convenient test naming
typedef ImplTest GetDataFromCheckpointTest;

TEST_F(GetDataFromCheckpointTest, ReturnsDefaultDataIfNoCpiFile)
{
    EXPECT_CALL(*input_, opt2bSet(_, _, _))
        .WillOnce(Return(false));
    DataFromCheckpoint actualData = impl_.getDataFromCheckpoint(true, false);
    EXPECT_EQ(defaultData_.partNumber_, actualData.partNumber_);
}

TEST_F(GetDataFromCheckpointTest, ReturnsDefaultIfNotSimMaster)
{
    EXPECT_CALL(*input_, opt2bSet(_, _, _))
        .WillOnce(Return(true));
    DataFromCheckpoint actualData = impl_.getDataFromCheckpoint(false, false);
    EXPECT_EQ(defaultData_.partNumber_, actualData.partNumber_);
}

TEST_F(GetDataFromCheckpointTest, ReturnsDefaultIfNeitherCheckpointOrLogFileExist)
{
    EXPECT_CALL(*input_, opt2bSet(_, _, _))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, ftp2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillOnce(Return(false))
        .WillOnce(Return(false));
    DataFromCheckpoint actualData = impl_.getDataFromCheckpoint(true, false);
    EXPECT_EQ(defaultData_.partNumber_, actualData.partNumber_);
}

TEST_F(GetDataFromCheckpointTest, ThrowsIfCheckpointFileDoesNotExistButLogFileDoes)
{
    EXPECT_CALL(*input_, opt2bSet(_, _, _))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, ftp2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillOnce(Return(false))
        .WillOnce(Return(true));
    EXPECT_THROW(impl_.getDataFromCheckpoint(true, false), InconsistentInputError);
}

TEST_F(GetDataFromCheckpointTest, ThrowsIfCheckpointFileCannotBeOpened)
{
    t_fileio *fp = NULL;

    EXPECT_CALL(*input_, opt2bSet(_, _, _))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _))
        .WillOnce(Return(fp));
    EXPECT_THROW(impl_.getDataFromCheckpoint(true, false), FileIOError);
}

TEST_F(GetDataFromCheckpointTest, ThrowsWhenCheckpointHadNoOutputFilenames)
{
    t_fileio             fp;
    int                  numFilenamesInCheckpoint = 0;
    gmx_file_position_t *empty                    = NULL;

    EXPECT_CALL(*input_, opt2bSet(_, _, _))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _))
        .WillOnce(Return(&fp));
    EXPECT_CALL(*input_, read_checkpoint_simulation_part_and_filenames(_, _, _, _))
        .WillOnce(DoAll(SetArgPointee<2>(numFilenamesInCheckpoint),
                        SetArgPointee<3>(empty)));
    EXPECT_THROW(impl_.getDataFromCheckpoint(true, false), InternalError);
}

TEST_F(GetDataFromCheckpointTest, ReadsCheckpointDataWhenAppropriate)
{
    t_fileio             fp;
    int                  checkpointPart = 10;
    gmx_file_position_t *empty          = NULL;

    EXPECT_CALL(*input_, opt2bSet(_, _, _))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _))
        .WillOnce(Return(&fp));
    EXPECT_CALL(*input_, read_checkpoint_simulation_part_and_filenames(_, _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(-1),
                        SetArgPointee<3>(empty)));
    DataFromCheckpoint actualData = impl_.getDataFromCheckpoint(true, false);
    EXPECT_EQ(checkpointPart, actualData.partNumber_);
}

// Set up data structures necessary for testing Impl::outputFilesExist
#ifndef DOXYGEN
const char                  *firstOutput[]  = { "firstOutput.log" };
const char                  *secondOutput[] = { "secondOutput.xvg" };
const char                  *input[]        = { "checkpoint.cpt" };
t_filenm                     fnm[]          = {
    { efLOG, "-o1",  "",  ffOPTWR, 1, const_cast<char **>(firstOutput) },
    { efXVG, "-o2",  "", ffOPTWR, 1, const_cast<char **>(secondOutput) },
    { efCPT, "-cpi", "", ffREAD | ffSET, 1, const_cast<char **>(input) },
};
ConstArrayRef<t_filenm>      fnms = fnm;
#endif      /* DOXYGEN */

//! Typedef for convenient test naming
typedef ImplTest OutputFileExistsTest;

TEST_F(OutputFileExistsTest, FindsFilesThatExist)
{
    impl_.mdrunFilenames_ = fnm;
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillRepeatedly(Return(true));
    EXPECT_TRUE(impl_.outputFileExists(firstOutput[0]));
    EXPECT_TRUE(impl_.outputFileExists(secondOutput[0]));
}

TEST_F(OutputFileExistsTest, CantFindFilesThatDontExist)
{
    impl_.mdrunFilenames_ = fnm;
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillRepeatedly(Return(false));
    EXPECT_FALSE(impl_.outputFileExists(secondOutput[0]));
}

TEST_F(OutputFileExistsTest, CantFindFilesThatAreNotInFnm)
{
    impl_.mdrunFilenames_ = fnm;
    EXPECT_FALSE(impl_.outputFileExists("missing"));
}

TEST_F(OutputFileExistsTest, CantFindInputFiles)
{
    impl_.mdrunFilenames_ = fnm;
    EXPECT_FALSE(impl_.outputFileExists(input[0]));
}

TEST_F(OutputFileExistsTest, CantFindFilesInEmptyList)
{
    impl_.mdrunFilenames_ = fnm;
    EXPECT_FALSE(impl_.outputFileExists(input[0]));
}

//! Test fixture for Impl::checkAllFilesForAppendingExist
class CheckAllFilesForAppendingExistTest : public ::testing::Test
{
    public:
        CheckAllFilesForAppendingExistTest() : impl_(NULL)
        {
            // The checkpoint file must name a log file
            expectedOutputFilenames_.push_back("md.log");
        }

        class MockImpl : public Impl
        {
            public:
                MockImpl(InputHandler *input) :
                    Impl(input, NULL, EmptyArrayRef())
                {
                }

                MOCK_CONST_METHOD1(outputFileExists, bool(const std::string&));
                MOCK_CONST_METHOD2(makeMessageWhenSetsOfFilesDontMatch,
                                   std::string(const std::vector<std::string> &,
                                               size_t));

        };

        //! Object under test
        MockImpl                 impl_;
        //! Input data for Impl::checkAllFilesForAppendingExist()
        std::vector<std::string> expectedOutputFilenames_;
};

TEST_F(CheckAllFilesForAppendingExistTest, ReturnsTrueWhenOutputFilesExist)
{
    EXPECT_CALL(impl_, outputFileExists(_))
        .WillRepeatedly(Return(true));
    EXPECT_NO_THROW(impl_.checkAllFilesForAppendingExist(expectedOutputFilenames_));
}

TEST_F(CheckAllFilesForAppendingExistTest, ReturnsFalseWhenNoOutputFilesExist)
{
    EXPECT_CALL(impl_, outputFileExists(_))
        .WillRepeatedly(Return(false));
    EXPECT_THROW(impl_.checkAllFilesForAppendingExist(expectedOutputFilenames_),
                 InconsistentInputError);
}

TEST_F(CheckAllFilesForAppendingExistTest, ThrowsWhenNotAllOutputFilesExist)
{
    // Add an output filename that will not exist
    expectedOutputFilenames_.push_back("ener.edr");
    EXPECT_CALL(impl_, outputFileExists(_))
        .WillOnce(Return(true))
        .WillRepeatedly(Return(false));
    EXPECT_CALL(impl_, makeMessageWhenSetsOfFilesDontMatch(_, _))
        .WillOnce(Return(""));
    EXPECT_THROW(impl_.checkAllFilesForAppendingExist(expectedOutputFilenames_),
                 InconsistentInputError);
}

//! Test fixture for Impl::canAppend
class CanAppendTest : public ::testing::Test
{
    public:
        CanAppendTest() : impl_(NULL)
        {
        }

        class MockImpl : public Impl
        {
            public:
                MockImpl(InputHandler *input) :
                    Impl(input, NULL, EmptyArrayRef())
                {
                }

                MOCK_CONST_METHOD1(checkAllFilesForAppendingExist,
                                   void(const std::vector<std::string> &expectedOutputFilenames));
        };

        //! Object under test
        MockImpl           impl_;
        //! Input data for Impl::canAppend()
        DataFromCheckpoint dataFromCheckpoint_;
};

TEST_F(CanAppendTest, AppendsWhenItShould)
{
    dataFromCheckpoint_.partNumber_ = 1;
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(_))
        .WillOnce(Return());
    dataFromCheckpoint_.expectedOutputFilenames_.push_back("md.log");
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
    EXPECT_FALSE(impl_.canAppend(dataFromCheckpoint_, false));
}

TEST_F(CanAppendTest, WontAppendIfCheckpointedSimulationUsedPartNumbers)
{
    dataFromCheckpoint_.partNumber_ = 1;
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(_))
        .WillOnce(Return());
    dataFromCheckpoint_.expectedOutputFilenames_.push_back("md.part0001.log");
    EXPECT_FALSE(impl_.canAppend(dataFromCheckpoint_, true));
}

TEST_F(CanAppendTest, ThrowsWhenCheckpointIsMalformed)
{
    dataFromCheckpoint_.partNumber_ = 1;
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(_))
        .WillOnce(Return());
    dataFromCheckpoint_.expectedOutputFilenames_.push_back("ener.edr");
    EXPECT_THROW(impl_.canAppend(dataFromCheckpoint_, true),
                 InternalError);
}

//! Test fixture for getRestartInformation
class GetRestartInformationTest : public ::testing::Test
{
    public:
        class MockImpl : public Impl
        {
            public:
                MockImpl(InputHandler *input, RankHandler *ranks) :
                    Impl(input, ranks, EmptyArrayRef())
                {
                }

                MOCK_CONST_METHOD1(checkAllFilesForAppendingExist,
                                   void(const std::vector<std::string> &expectedOutputFilenames));
        };

        GetRestartInformationTest() : input_(new TestInputHandler),
                                      ranks_(new TestRankHandler),
                                      impl_(input_, ranks_),
                                      fp_(),
                                      numFilenamesInCheckpoint_(1),
                                      outputFilenames_()
        {
            snew(outputFilenames_, numFilenamesInCheckpoint_);
            strncpy(outputFilenames_[0].filename, "md.log", STRLEN);
        }

        /*! \brief Handle to mock input object; impl_ manages its
            lifetime */
        TestInputHandler    *input_;
        /*! \brief Handle to mock object for handling multiple ranks;
            impl_ manages its lifetime */
        TestRankHandler     *ranks_;
        //! Object under test
        MockImpl             impl_;
        //! Fake file object opened by gmx_fio_open()
        t_fileio             fp_;
        //! Output of read_checkpoint_simulation_part_and_filenames
        int                  numFilenamesInCheckpoint_;
        //! Output of read_checkpoint_simulation_part_and_filenames
        gmx_file_position_t *outputFilenames_;
};

TEST_F(GetRestartInformationTest, DoesPlainStartFromCptIfNoCpiOptionUsed)
{
    // Mock behaviour with multiple ranks
    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(true));
    EXPECT_CALL(*ranks_, hasMultipleRanks()).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, checkAcrossMultiSim(_, _, _, _)).WillOnce(Return());

    // Mock checking for checkpoint file option
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(false));

    // Call the code to test
    RestartInformation info = impl_.getRestartInformation(true);

    // Assert postconditions
    EXPECT_FALSE(info.bWillStartFromCpt_);
    EXPECT_FALSE(info.bWillAppendFiles_);
}

TEST_F(GetRestartInformationTest, DoesPlainStartIfNeitherCptOrLogFileExists)
{
    // Mock behaviour with multiple ranks
    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    // The next use of false just suppresses output to stdout
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, hasMultipleRanks()).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, checkAcrossMultiSim(_, _, _, _)).WillOnce(Return());

    // Mock checking for checkpoint and log files
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillOnce(Return(false))
        .WillOnce(Return(false));
    EXPECT_CALL(*input_, ftp2fn(_, _, _));

    // Call the code to test
    RestartInformation info = impl_.getRestartInformation(true);

    // Assert postconditions
    EXPECT_FALSE(info.bWillStartFromCpt_);
    EXPECT_FALSE(info.bWillAppendFiles_);
}

TEST_F(GetRestartInformationTest, ThrowsIfCptDoesNotExistAndLogFileDoesExist)
{
    // Mock behaviour with multiple ranks
    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(true));
    // Mock checking for checkpoint and log files
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_))
        .WillOnce(Return(false))
        .WillOnce(Return(true));
    EXPECT_CALL(*input_, ftp2fn(_, _, _));

    // Call the code to test
    EXPECT_THROW(impl_.getRestartInformation(true), InconsistentInputError);
}

TEST_F(GetRestartInformationTest, ThrowsIfCptCannotBeOpened)
{
    t_fileio *nullPointer = NULL;

    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(true));
    // Mock checking for checkpoint file
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_)).WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _)).WillOnce(Return(nullPointer));

    // Call the code to test
    EXPECT_THROW(impl_.getRestartInformation(true), FileIOError);
}

TEST_F(GetRestartInformationTest, ThrowsIfCptFailsConsistencyChecks)
{
    int checkpointPart = 10;

    // Mock behaviour with multiple ranks
    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(true));
    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_)).WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _)).WillOnce(Return(&fp_));
    EXPECT_CALL(*input_, read_checkpoint_simulation_part_and_filenames(_, _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));
    // Mock consistency checks
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(_))
        .WillOnce(Throw(InconsistentInputError("")));

    // Call the code to test
    EXPECT_THROW(impl_.getRestartInformation(true), InconsistentInputError);
}

TEST_F(GetRestartInformationTest, StartsFromCptAndAppendIfEverythingIsOK)
{
    int checkpointPart = 10;

    // Mock behaviour with multiple ranks
    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(true));
    EXPECT_CALL(*ranks_, hasMultipleRanks()).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, checkAcrossMultiSim(_, _, _, _)).WillOnce(Return());

    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_)).WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _)).WillOnce(Return(&fp_));
    EXPECT_CALL(*input_, read_checkpoint_simulation_part_and_filenames(_, _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));
    // Mock consistency checks
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(_)).WillOnce(Return());

    // Call the code to test
    RestartInformation info = impl_.getRestartInformation(true);

    // Assert postconditions
    EXPECT_TRUE(info.bWillStartFromCpt_);
    EXPECT_TRUE(info.bWillAppendFiles_);
}

TEST_F(GetRestartInformationTest, StartsFromCptAndUsesPartNumbersIfOldRunDid)
{
    int checkpointPart = 10;
    strncpy(outputFilenames_[0].filename, "md.part0010.log", STRLEN);

    // Mock behaviour with multiple ranks
    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    // The next use of false just suppresses output to stdout
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(true)).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, hasMultipleRanks()).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, checkAcrossMultiSim(_, _, _, _)).WillOnce(Return());

    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_)).WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _)).WillOnce(Return(&fp_));
    EXPECT_CALL(*input_, read_checkpoint_simulation_part_and_filenames(_, _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));
    // Mock consistency checks
    EXPECT_CALL(impl_, checkAllFilesForAppendingExist(_)).WillOnce(Return());

    // Call the code to test
    RestartInformation info = impl_.getRestartInformation(true);

    // Assert postconditions
    EXPECT_TRUE(info.bWillStartFromCpt_);
    EXPECT_FALSE(info.bWillAppendFiles_);
}

TEST_F(GetRestartInformationTest, StartsFromCptAndUsesPartNumbersIfNoAppendIsSet)
{
    int checkpointPart = 10;

    // Mock behaviour with multiple ranks
    EXPECT_CALL(*ranks_, isSimMaster()).WillOnce(Return(true));
    // The next use of false just suppresses output to stdout
    EXPECT_CALL(*ranks_, isMultiMaster()).WillOnce(Return(true)).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, hasMultipleRanks()).WillOnce(Return(false));
    EXPECT_CALL(*ranks_, checkAcrossMultiSim(_, _, _, _)).WillOnce(Return());

    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(*input_, opt2bSet(_, _, _)).WillOnce(Return(true));
    EXPECT_CALL(*input_, opt2fn(_, _, _));
    EXPECT_CALL(*input_, gmx_fexist(_)).WillOnce(Return(true));
    EXPECT_CALL(*input_, gmx_fio_open(_, _)).WillOnce(Return(&fp_));
    EXPECT_CALL(*input_, read_checkpoint_simulation_part_and_filenames(_, _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));

    // Call the code to test
    RestartInformation info = impl_.getRestartInformation(false);

    // Assert postconditions
    EXPECT_TRUE(info.bWillStartFromCpt_);
    EXPECT_FALSE(info.bWillAppendFiles_);
}

} // namespace

} // namespace

} // namespace
