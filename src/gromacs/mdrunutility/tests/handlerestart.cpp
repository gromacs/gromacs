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
 * Implements tests for handling mdrun restarts
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/mdrunutility/handlerestart-impl.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/mpihandlerinterface.h"
#include "gromacs/utility/smalloc.h"

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

class MockMpiHandler : public MpiHandlerInterface
{
    public:
        MockMpiHandler() {};

        MOCK_CONST_METHOD0(hasMultipleRanks, bool());
        MOCK_CONST_METHOD0(isMaster, bool());
        MOCK_CONST_METHOD0(isSimMaster, bool());
        MOCK_CONST_METHOD0(isMultiSim, bool());
        MOCK_CONST_METHOD0(isMultiMaster, bool());
        MOCK_CONST_METHOD2(broadcast, void(int, void *));
        //! Empty implementation because we always want to skip this check
        virtual void checkAcrossMultiSim(FILE *, int, const char *, bool) const {};
};

/* Ideally, these would be static const members of
   HandleRestartTest, but C++98 can't cope with in-class array
   initializers, which means that the ArrayRef template constructors
   that are based on these arrays can't work. */
//! Fake mdrun argv array, used in most tests
const char *HandleRestartTestArgv_g[] = { "mdrun", "-cpi", "state.cpt", "-g", "md.log" };
//! Fake mdrun filenm structure to fill from argv, used in most tests
t_filenm    HandleRestartTestFilenames_g[] = { { efCPT, "-cpi", NULL, ffWRITE },
                                               { efLOG, "-g", NULL, ffWRITE } };

//! Test fixture for getRestartInformation
class HandleRestartTest : public ::testing::Test
{
    public:
        HandleRestartTest() : fileIOHandler_(),
                              mpiHandler_(),
                              numFilenamesInCheckpoint_(1),
                              outputFilenames_()
        {
            snew(outputFilenames_, numFilenamesInCheckpoint_);
            strncpy(outputFilenames_[0].filename, "md.log", STRLEN);
        }
        ~HandleRestartTest()
        {
            // TODO We ought to need to free outputFilenames_, but
            // doing so causes a memory error in apparently unrelated
            // code
            //sfree(outputFilenames_);
        }
        //! Called once for all tests in test case
        static void SetUpTestCase()
        {
            /* Build the rest of the filenames array from the argv
               inputs. Only needs to be done once for all test cases,
               because the data is logically const. */
            parseCommonArgs(HandleRestartTestArgv_g,
                            HandleRestartTestFilenames_g);
        }

        //! Mock input object
        MockFileIOHandler    fileIOHandler_;
        //! Mock object for handling multiple ranks
        MockMpiHandler       mpiHandler_;
        //! Number of names of output files recorded in the checkpoint
        int                  numFilenamesInCheckpoint_;
        //! Array of structures containing names of output files recorded in the checkpoint
        gmx_file_position_t *outputFilenames_;
};

TEST_F(HandleRestartTest, DoesPlainStartFromCptIfNoCpiOptionWasGiven)
{
    // Mock behaviour with multiple ranks
    EXPECT_CALL(mpiHandler_, isSimMaster()).WillOnce(Return(true));
    // The next use of false just suppresses output to stdout
    EXPECT_CALL(mpiHandler_, isMultiMaster()).WillOnce(Return(false));
    EXPECT_CALL(mpiHandler_, hasMultipleRanks()).WillOnce(Return(false));

    /* Call the code to test, but use an empty array for filenames so
     * that no -cpi is given */
    RestartInformation info =
        gmx::handleRestart(&fileIOHandler_, &mpiHandler_, EmptyArrayRef(), true);

    // Assert postconditions
    EXPECT_FALSE(info.bStartFromCpt_);
    EXPECT_FALSE(info.bAppendFiles_);
}

TEST_F(HandleRestartTest, DoesPlainStartIfNeitherCptOrLogFileExists)
{
    // Mock behaviour with multiple ranks
    EXPECT_CALL(mpiHandler_, isSimMaster()).WillOnce(Return(true));
    // The next use of false just suppresses output to stdout
    EXPECT_CALL(mpiHandler_, isMultiMaster()).WillOnce(Return(false));
    EXPECT_CALL(mpiHandler_, hasMultipleRanks()).WillOnce(Return(false));

    // Mock checking for checkpoint and log files
    EXPECT_CALL(fileIOHandler_, gmx_fexist(_))
        .WillOnce(Return(false))
        .WillOnce(Return(false));

    // Call the code to test, using the normal filenames
    RestartInformation info =
        gmx::handleRestart(&fileIOHandler_, &mpiHandler_,
                           CommandLineFilenames(HandleRestartTestFilenames_g), true);

    // Assert postconditions
    EXPECT_FALSE(info.bStartFromCpt_);
    EXPECT_FALSE(info.bAppendFiles_);
}

TEST_F(HandleRestartTest, ThrowsIfCptDoesNotExistAndLogFileDoesExist)
{
    // Mock behaviour with multiple ranks
    EXPECT_CALL(mpiHandler_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(mpiHandler_, isMultiMaster()).WillOnce(Return(true));
    // Mock checking for checkpoint and log files
    EXPECT_CALL(fileIOHandler_, gmx_fexist(_))
        .WillOnce(Return(false))
        .WillOnce(Return(true));

    // Call the code to test, using the normal filenames
    EXPECT_THROW_GMX(gmx::handleRestart(&fileIOHandler_, &mpiHandler_,
                                        CommandLineFilenames(HandleRestartTestFilenames_g), true), InconsistentInputError);
}

TEST_F(HandleRestartTest, ThrowsIfCptFailsConsistencyChecks)
{
    int checkpointPart = 10;

    /* Add another output filename to the set from the checkpoint
       file, so that we get inconsistency with the command line, and
       the exception will be thrown. */
    srenew(outputFilenames_, numFilenamesInCheckpoint_+1);
    strncpy(outputFilenames_[numFilenamesInCheckpoint_].filename, "dummy.xvg", STRLEN);
    numFilenamesInCheckpoint_++;

    // Mock behaviour with multiple ranks
    EXPECT_CALL(mpiHandler_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(mpiHandler_, isMultiMaster()).WillOnce(Return(true));
    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(fileIOHandler_, gmx_fexist(StrEq("md.log"))).Times(3).WillRepeatedly(Return(true));
    EXPECT_CALL(fileIOHandler_, gmx_fexist(StrEq("state.cpt"))).Times(1).WillRepeatedly(Return(true));
    EXPECT_CALL(fileIOHandler_, read_checkpoint_simulation_part_and_filenames(StrEq("state.cpt"), _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));

    // Call the code to test, using the normal filenames in impl_
    EXPECT_THROW_GMX(gmx::handleRestart(&fileIOHandler_, &mpiHandler_,
                                        CommandLineFilenames(HandleRestartTestFilenames_g), true), InconsistentInputError);
}

TEST_F(HandleRestartTest, StartsFromCptAndAppendIfEverythingIsOK)
{
    int checkpointPart = 10;

    // Mock behaviour with multiple ranks
    EXPECT_CALL(mpiHandler_, isSimMaster()).WillOnce(Return(true));
    // These use of false just suppress output to stdout
    EXPECT_CALL(mpiHandler_, isMultiMaster()).Times(2).WillRepeatedly(Return(false));
    EXPECT_CALL(mpiHandler_, hasMultipleRanks()).WillOnce(Return(false));

    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(fileIOHandler_, gmx_fexist(StrEq("md.log"))).Times(1).WillRepeatedly(Return(true));
    EXPECT_CALL(fileIOHandler_, gmx_fexist(StrEq("state.cpt"))).Times(1).WillRepeatedly(Return(true));
    EXPECT_CALL(fileIOHandler_, read_checkpoint_simulation_part_and_filenames(StrEq("state.cpt"), _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));

    // Call the code to test, using the normal filenames
    RestartInformation info =
        gmx::handleRestart(&fileIOHandler_, &mpiHandler_,
                           CommandLineFilenames(HandleRestartTestFilenames_g), true);

    // Assert postconditions
    EXPECT_TRUE(info.bStartFromCpt_);
    EXPECT_TRUE(info.bAppendFiles_);
}

TEST_F(HandleRestartTest, StartsFromCptAndUsesPartNumbersIfOldRunDid)
{
    int checkpointPart = 10;
    // Fake mdrun inputs for this test need to follow the *.part%04d.* pattern
    strncpy(outputFilenames_[0].filename, "md.part0010.log", STRLEN);
    const char *argv[]      = { "mdrun", "-cpi", "state.cpt", "-g", "md.part0010.log" };
    t_filenm    filenames[] = { { efCPT, "-cpi", NULL, ffWRITE },
                                { efLOG, "-g", NULL, ffWRITE } };
    // parse that argv into the filename data structure
    parseCommonArgs(argv, filenames);

    // Mock behaviour with multiple ranks
    EXPECT_CALL(mpiHandler_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(mpiHandler_, isMultiMaster()).WillOnce(Return(true));
    EXPECT_CALL(mpiHandler_, hasMultipleRanks()).WillOnce(Return(false));

    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(fileIOHandler_, gmx_fexist(StrEq("state.cpt"))).WillOnce(Return(true));
    EXPECT_CALL(fileIOHandler_, gmx_fexist(StrEq("md.part0010.log"))).WillOnce(Return(true));
    EXPECT_CALL(fileIOHandler_, read_checkpoint_simulation_part_and_filenames(StrEq("state.cpt"), _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));

    // Call the code to test
    RestartInformation info =
        gmx::handleRestart(&fileIOHandler_, &mpiHandler_,
                           CommandLineFilenames(filenames), true);

    // Assert postconditions
    EXPECT_TRUE(info.bStartFromCpt_);
    EXPECT_FALSE(info.bAppendFiles_);
}

TEST_F(HandleRestartTest, StartsFromCptAndUsesPartNumbersIfNoAppendIsSet)
{
    int checkpointPart = 10;

    // Mock behaviour with multiple ranks
    EXPECT_CALL(mpiHandler_, isSimMaster()).WillOnce(Return(true));
    EXPECT_CALL(mpiHandler_, isMultiMaster()).WillOnce(Return(true));
    EXPECT_CALL(mpiHandler_, hasMultipleRanks()).WillOnce(Return(false));

    // Mock checking for and reading of checkpoint file
    EXPECT_CALL(fileIOHandler_, gmx_fexist(_)).WillOnce(Return(true));
    EXPECT_CALL(fileIOHandler_, read_checkpoint_simulation_part_and_filenames(StrEq("state.cpt"), _, _, _))
        .WillOnce(DoAll(SetArgPointee<1>(checkpointPart),
                        SetArgPointee<2>(numFilenamesInCheckpoint_),
                        SetArgPointee<3>(outputFilenames_)));

    // Call the code to test, using the normal filenames
    RestartInformation info =
        gmx::handleRestart(&fileIOHandler_, &mpiHandler_,
                           CommandLineFilenames(HandleRestartTestFilenames_g), false);

    // Assert postconditions
    EXPECT_TRUE(info.bStartFromCpt_);
    EXPECT_FALSE(info.bAppendFiles_);
}

} // namespace

} // namespace

} // namespace
