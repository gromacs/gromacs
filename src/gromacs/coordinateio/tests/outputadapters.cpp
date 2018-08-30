/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!\file
 * \internal
 * \brief
 * Tests for outputmanager
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */

#include "gmxpre.h"

#include "outputadapters.h"

namespace gmx
{

namespace test
{

TEST_P(SetAtomsSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(SetAtomsUnSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(AnyOutputSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_F(OutputSelectorDeathTest, RejectsBadSelection)
{
    prepareTest();
}

TEST_P(SetVelocitySupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(SetVelocityUnSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(SetForceSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(SetForceUnSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(SetPrecisionSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(SetPrecisionUnSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST_P(OutputTurnedOffSupportedFiles, Works)
{
    prepareTest(GetParam());
}

TEST(SetTime, UnchangedTimeWorks)
{
    SetTime    testTime(0, 1, ChangeFrameTimeType::efUnchanged);

    t_trxframe testFrame;
    clear_trxframe(&testFrame, true);
    testFrame.time = 23;

    testTime.processFrame(1, &testFrame);

    EXPECT_EQ(23, testFrame.time);

    testTime.processFrame(2, &testFrame);

    EXPECT_EQ(23, testFrame.time);
}

TEST(SetTime, ChangeTimeStepWorks)
{
    SetTime    testTime(0, 1, ChangeFrameTimeType::efTimeStep);

    t_trxframe testFrame;
    clear_trxframe(&testFrame, true);
    testFrame.time = 23;

    testTime.processFrame(0, &testFrame);

    EXPECT_EQ(23, testFrame.time);

    testTime.processFrame(1, &testFrame);

    EXPECT_EQ(24, testFrame.time);
}

TEST(SetTime, ChangeStartTimeWorks)
{
    SetTime    testTime(0, 1, ChangeFrameTimeType::efStartTime);

    t_trxframe testFrame;
    clear_trxframe(&testFrame, true);
    testFrame.time = 23;

    testTime.processFrame(0, &testFrame);

    EXPECT_EQ(0, testFrame.time);

    testFrame.time = 24;

    testTime.processFrame(1, &testFrame);

    EXPECT_EQ(1, testFrame.time);
}

TEST(SetTime, ChangeBothWorks)
{
    SetTime    testTime(0, 1, ChangeFrameTimeType::efBothTime);

    t_trxframe testFrame;
    clear_trxframe(&testFrame, true);
    testFrame.time = 23;

    testTime.processFrame(0, &testFrame);

    EXPECT_EQ(0, testFrame.time);

    testFrame.time = 26;

    testTime.processFrame(1, &testFrame);

    EXPECT_EQ(1, testFrame.time);
}


INSTANTIATE_TEST_CASE_P(ModuleSupported,
                        SetAtomsSupportedFiles, ::testing::ValuesIn(setAtomsSupported));

INSTANTIATE_TEST_CASE_P(ModuleUnSupported,
                        SetAtomsUnSupportedFiles, ::testing::ValuesIn(setAtomsUnSupported));

INSTANTIATE_TEST_CASE_P(ModuleSupported,
                        AnyOutputSupportedFiles, ::testing::ValuesIn(anySupported));

INSTANTIATE_TEST_CASE_P(ModuleSupported,
                        SetVelocitySupportedFiles, ::testing::ValuesIn(setVelocitySupported));

INSTANTIATE_TEST_CASE_P(ModuleUnSupported,
                        SetVelocityUnSupportedFiles, ::testing::ValuesIn(setVelocityUnSupported));

INSTANTIATE_TEST_CASE_P(ModuleSupported,
                        SetForceSupportedFiles, ::testing::ValuesIn(setForceSupported));

INSTANTIATE_TEST_CASE_P(ModuleUnSupported,
                        SetForceUnSupportedFiles, ::testing::ValuesIn(setForceUnSupported));

INSTANTIATE_TEST_CASE_P(ModuleSupported,
                        SetPrecisionSupportedFiles, ::testing::ValuesIn(setPrecisionSupported));

INSTANTIATE_TEST_CASE_P(ModuleUnSupported,
                        SetPrecisionUnSupportedFiles, ::testing::ValuesIn(setPrecisionUnSupported));

INSTANTIATE_TEST_CASE_P(ModuleSupported,
                        OutputTurnedOffSupportedFiles, ::testing::ValuesIn(anySupported));
} // namespace test

} // namespace gmx
