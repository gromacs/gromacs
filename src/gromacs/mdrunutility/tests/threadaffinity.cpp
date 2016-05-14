/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/mdrunutility/threadaffinity.h"

#include <memory>

#include <gmock/gmock.h>

#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/smalloc.h"

#include "threadaffinitytest.h"

namespace
{

class ThreadAffinityTest : public ::testing::Test
{
    public:
        ThreadAffinityTest()
        {
        }
        ~ThreadAffinityTest()
        {
            sfree(cr_);
            sfree(hwOpt_);
        }

        void setTestParameters(int affinityOption, int logicalProcessorCount)
        {
            snew(cr_, 1);
            cr_->nnodes = 1;
            cr_->nodeid = 0;
            snew(hwOpt_, 1);
            hwOpt_->thread_affinity = affinityOption;
            hwTop_.reset(new gmx::HardwareTopology(logicalProcessorCount));
        }

        void setAffinity(int nthread_local)
        {
            gmx_set_thread_affinity(nullptr, cr_, hwOpt_, *hwTop_,
                                    nthread_local, &affinityAccess_);
        }

        gmx::test::MockThreadAffinityAccess     affinityAccess_;
        gmx_hw_opt_t                           *hwOpt_;

    private:
        t_commrec                              *cr_;
        std::unique_ptr<gmx::HardwareTopology>  hwTop_;
};

TEST_F(ThreadAffinityTest, DoesNothingWhenDisabled)
{
    setTestParameters(threadaffOFF, 1);
    setAffinity(1);
}

TEST_F(ThreadAffinityTest, DoesNothingWhenNotSupported)
{
    using ::testing::Return;
    setTestParameters(threadaffAUTO, 1);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported())
        .WillOnce(Return(false));
    setAffinity(1);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooFewThreads)
{
    setTestParameters(threadaffAUTO, 4);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    setAffinity(2);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooManyThreads)
{
    setTestParameters(threadaffAUTO, 4);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithUnknownHardware)
{
    setTestParameters(threadaffON, 0);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    setAffinity(2);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooManyThreads)
{
    setTestParameters(threadaffON, 4);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooLargeOffset)
{
    setTestParameters(threadaffON, 4);
    hwOpt_->core_pinning_offset = 2;
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    setAffinity(3);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooLargeStride)
{
    setTestParameters(threadaffON, 4);
    hwOpt_->core_pinning_stride = 2;
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    setAffinity(3);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWithAuto)
{
    setTestParameters(threadaffAUTO, 1);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(0));
    setAffinity(1);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWhenForced)
{
    setTestParameters(threadaffON, 2);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(0));
    setAffinity(1);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWithOffsetWhenForced)
{
    setTestParameters(threadaffON, 4);
    hwOpt_->core_pinning_offset = 2;
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(2));
    setAffinity(1);
}

TEST_F(ThreadAffinityTest, HandlesPinningFailureWithSingleThread)
{
    using ::testing::Return;
    setTestParameters(threadaffAUTO, 1);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(0))
        .WillOnce(Return(false));
    setAffinity(1);
}

// TODO: If it wouldn't result in a multitude of #if's, it would be nice
// to somehow indicate in a no-OpenMP build that some tests are missing.
#if GMX_OPENMP
TEST_F(ThreadAffinityTest, PinsMultipleThreadsWithAuto)
{
    setTestParameters(threadaffAUTO, 2);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(0));
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(1));
    setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsMultipleThreadsWithStrideWhenForced)
{
    setTestParameters(threadaffON, 4);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    hwOpt_->core_pinning_stride = 2;
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(0));
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(2));
    setAffinity(2);
}

TEST_F(ThreadAffinityTest, HandlesPinningFailureWithOneThreadFailing)
{
    using ::testing::Return;
    setTestParameters(threadaffON, 2);
    EXPECT_CALL(affinityAccess_, isThreadAffinitySupported());
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(0));
    EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(1))
        .WillOnce(Return(false));
    setAffinity(2);
}
#endif

} // namespace
