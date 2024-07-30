/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/mdrunutility/threadaffinity.h"

#include "config.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/hardware/hw_info.h"

#include "threadaffinitytest.h"

namespace gmx
{
namespace test
{
namespace
{

class ThreadAffinityTest : public ::testing::Test
{
public:
    gmx::test::ThreadAffinityTestHelper helper_;
};

TEST_F(ThreadAffinityTest, DoesNothingWhenDisabled)
{
    helper_.setAffinityOption(ThreadAffinity::Off);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, DoesNothingWhenNotSupported)
{
    helper_.setAffinitySupported(false);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooFewUserSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("The number of threads is not equal to the number of");
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooManyUserSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper_.setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooManyAutoSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.setTotNumThreadsIsAuto(true);
    helper_.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper_.setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithUnknownHardware)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(0);
    helper_.expectWarningMatchingRegex("No information on available logical cpus");
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooManyThreads)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper_.setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooLargeOffset)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(2, 0);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Applying core pinning offset 2");
    helper_.expectWarningMatchingRegex("Requested offset too large");
    helper_.setAffinity(3);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooLargeStride)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(0, 2);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Requested stride too large");
    helper_.setAffinity(3);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWithAuto)
{
    helper_.setLogicalProcessorCount(1);
    helper_.expectAffinitySet(0);
    helper_.expectPinningMessage(false, 1);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWhenForced)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(2);
    helper_.expectPinningMessage(false, 2);
    helper_.expectAffinitySet(0);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWithOffsetWhenForced)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(2, 0);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Applying core pinning offset 2");
    helper_.expectPinningMessage(false, 2);
    helper_.expectAffinitySet(2);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, HandlesPinningFailureWithSingleThread)
{
    helper_.setLogicalProcessorCount(1);
    helper_.expectPinningMessage(false, 1);
    // It would be good to check for the warning, but that currently goes to stderr
    helper_.expectAffinitySetThatFails(0);
    helper_.setAffinity(1);
}

// TODO: If it wouldn't result in a multitude of #if's, it would be nice
// to somehow indicate in a no-OpenMP build that some tests are missing.
#if GMX_OPENMP
TEST_F(ThreadAffinityTest, PinsMultipleThreadsWithAuto)
{
    helper_.setLogicalProcessorCount(2);
    helper_.expectPinningMessage(false, 1);
    helper_.expectAffinitySet({ 0, 1 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsMultipleThreadsWithStrideWhenForced)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(0, 2);
    helper_.setLogicalProcessorCount(4);
    helper_.expectPinningMessage(true, 2);
    helper_.expectAffinitySet({ 0, 2 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsWithAutoAndFewerAutoSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.setTotNumThreadsIsAuto(true);
    helper_.expectPinningMessage(false, 2);
    helper_.expectAffinitySet({ 0, 2 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, HandlesPinningFailureWithOneThreadFailing)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(2);
    helper_.expectPinningMessage(false, 1);
    helper_.expectGenericFailureMessage();
    helper_.expectAffinitySet(0);
    helper_.expectAffinitySetThatFails(1);
    helper_.setAffinity(2);
}
#endif

} // namespace
} // namespace test
} // namespace gmx
