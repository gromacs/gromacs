/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Tests for gmx::Mutex
 *
 * These tests ensure that basic mutual-exclusion properties hold.
 * Note that no testing can prove there isn't a bug, but if one
 * exists, then these tests might expose one.
 *
 * In particular, try_lock can be implemented differently on different
 * platforms, or with different default mutex types, so we should
 * check that the behaviour continues to conform with the thread-MPI
 * documentation.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "gromacs/utility/mutex.h"

#include "config.h"

#include <future>
#include <gtest/gtest.h>


namespace gmx
{
namespace test
{
namespace
{

//! Convenience definition.
using Lock = gmx::lock_guard<Mutex>;

TEST(MutexBasicTest, CanBeMade)
{
    Mutex m;
}

TEST(MutexBasicTest, CanBeLocked)
{
    Mutex m;
    ASSERT_NO_THROW(m.lock());
    m.unlock();
}

TEST(MutexBasicTest, CanBeTryLocked)
{
    Mutex m;
    ASSERT_TRUE(m.try_lock());
    m.unlock();
}

TEST(MutexBasicTest, CanBeUsedInLockGuard)
{
    Mutex m;
    Lock  g(m);
}

//! A shared value for a mutex to protect
static int   g_sharedValue;
//! A mutex to protect a shared value
static Mutex g_sharedValueMutex;

//! Function type for asynchronous tasks.
using TaskType = std::function<int(void)>;

//! A task that just does work.
int updateSharedValue()
{
    return ++g_sharedValue;
}

//! A task that does work after it gets the mutex.
int updateSharedValueWithLock()
{
    Lock guard(g_sharedValueMutex);
    return updateSharedValue();
}

//! A task that does work only if it can get the mutex immediately.
int updateSharedValueWithTryLock()
{
    // Special return value to signal when work was not done because
    // the lock was not acquired.
    int result = -1;
    if (g_sharedValueMutex.try_lock())
    {
        result = updateSharedValue();
        g_sharedValueMutex.unlock();
    }
    return result;
}

/*! \brief Parameterized test fixture.
 *
 * Checks that different launch policies work. In further tests of
 * mutual exclusion, we need to specify std::thread::async, to require
 * that a thread actually launched. The default policy permits the
 * std:: implementation to avoid launching a thread, and at least the
 * behaviour of thread-MPI try_lock also varies with the threading
 * implementation underlying it. */
class DifferentTasksTest : public ::testing::TestWithParam<TaskType>
{
    public:
        //! Constructor
        DifferentTasksTest()
        {
            g_sharedValue = 0;
        }
        //! Check the results
        void checkResults()
        {
            int result;
            EXPECT_NO_THROW(result = futureResult_.get()) << "Future should not contain an exception";
            EXPECT_EQ(1, result) << "Task should have run";
            EXPECT_EQ(1, g_sharedValue) << "Shared value should be updated";
        }
        //! Contains the result the task returns.
        std::future<int> futureResult_;
};

TEST_P(DifferentTasksTest, StdAsyncWorksWithDefaultPolicy)
{
    auto task = GetParam();
    EXPECT_NO_THROW(futureResult_ = std::async(task)) << "Async should succeed";
    checkResults();
}

TEST_P(DifferentTasksTest, StdAsyncWorksWithAsyncLaunchPolicy)
{
    auto task = GetParam();
    EXPECT_NO_THROW(futureResult_ = std::async(std::launch::async, task)) << "Async should succeed";
    checkResults();
}

TEST_P(DifferentTasksTest, StdAsyncWorksWithDeferredLaunchPolicy)
{
    auto task = GetParam();
    EXPECT_NO_THROW(futureResult_ = std::async(std::launch::deferred, task)) << "Async should succeed";
    checkResults();
}

// Test that the different launch policies work with the different tasks
INSTANTIATE_TEST_CASE_P(WithAndWithoutMutex, DifferentTasksTest, ::testing::Values(updateSharedValue, updateSharedValueWithLock, updateSharedValueWithTryLock));

TEST(MutexTaskTest, MutualExclusionWorksWithLock)
{
    g_sharedValue = 0;
    std::future<int> result;
    {
        // Hold the mutex, launch a lock attempt on another
        // thread, check that the shared value isn't changed, then
        // release the mutex by leaving the scope, after which the
        // other thread's lock can get the mutex.
        Lock guard(g_sharedValueMutex);
        result = std::async(std::launch::async, updateSharedValueWithLock);
        EXPECT_EQ(0, g_sharedValue) << "Task should not have run yet";
    }
    EXPECT_EQ(1, result.get()) << "Task should have run";
    EXPECT_EQ(1, g_sharedValue) << "Shared value should be updated";
}

TEST(MutexTaskTest, MutualExclusionWorksWithTryLockOnOtherThread)
{
    g_sharedValue = 0;
    {
        // Hold the mutex, launch a try_lock attempt on another
        // thread, check that the shared value isn't changed, then
        // make sure the try_lock attempt has returned, double check
        // that the shared value isn't changed, and release the mutex
        // by leaving the scope.
        Lock guard(g_sharedValueMutex);
        auto result = std::async(std::launch::async, updateSharedValueWithTryLock);
        EXPECT_EQ(0, g_sharedValue) << "Data race detected";
        EXPECT_EQ(-1, result.get()) << "The try_lock should fail";
        EXPECT_EQ(0, g_sharedValue) << "Task should not have run";
    }
    EXPECT_EQ(0, g_sharedValue) << "Mutex release can't affect the protected value";
}

TEST(MutexTaskTest, MutualExclusionWorksWithTryLockOnSameThread)
{
    g_sharedValue = 0;
    int finalSharedValue = GMX_NATIVE_WINDOWS ? 1 : 0;
    {
        // Hold the mutex and launch a try_lock attempt on this
        // thread. Behaviour then varies with the implementation
        // underlying thread-MPI.
        Lock guard(g_sharedValueMutex);
        int  result = updateSharedValueWithTryLock();
        if (GMX_NATIVE_WINDOWS)
        {
            EXPECT_EQ(1, result) << "The try_lock should succeed";
            EXPECT_EQ(finalSharedValue, g_sharedValue) << "Task should have run";
        }
        else
        {
            EXPECT_EQ(-1, result) << "The try_lock should fail";
            EXPECT_EQ(finalSharedValue, g_sharedValue) << "Task should not have run";
        }
    }
    EXPECT_EQ(finalSharedValue, g_sharedValue) << "Mutex release can't affect the protected value";
}

} // namespace
} // namespace
} // namespace
