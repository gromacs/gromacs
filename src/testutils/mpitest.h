/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Helper functions for MPI tests to make thread-MPI look like real MPI.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_MPITEST_H
#define GMX_TESTUTILS_MPITEST_H

#include "config.h"

#include <functional>
#include <type_traits>

#include "gromacs/utility/basenetwork.h"

namespace gmx
{
namespace test
{

/*! \brief
 * Returns the number of MPI ranks to use for an MPI test.
 *
 * For thread-MPI builds, this will return the requested number of ranks
 * even before the thread-MPI threads have been started.
 *
 * \ingroup module_testutils
 */
int getNumberOfTestMpiRanks();
//! \cond internal
/*! \brief
 * Helper function for GMX_MPI_TEST().
 *
 * \ingroup module_testutils
 */
bool threadMpiTestRunner(std::function<void()> testBody);
//! \endcond

/*! \brief
 * Declares that this test is an MPI-enabled unit test.
 *
 * \param[in] expectedRankCount Expected number of ranks for this test.
 *     The test will fail if run with unsupported number of ranks.
 *
 * To write unit tests that run under MPI, you need to do a few things:
 *  - Put GMX_MPI_TEST() as the first statement in your test body and
 *    specify the number of ranks this test expects.
 *  - Declare your unit test in CMake with gmx_add_mpi_unit_test().
 *    Note that all tests in the binary should fulfill the conditions above,
 *    and work with the same number of ranks.
 * TODO: Figure out a mechanism for mixing tests with different rank counts in
 * the same binary (possibly, also MPI and non-MPI tests).
 *
 * When you do the above, the following will happen:
 *  - The test will get compiled only if thread-MPI or real MPI is enabled.
 *  - The test will get executed on the number of ranks specified.
 *    If you are using real MPI, the whole test binary is run under MPI and
 *    test execution across the processes is synchronized (GMX_MPI_TEST()
 *    actually has no effect in this case, the synchronization is handled at a
 *    higher level).
 *    If you are using thread-MPI, GMX_MPI_TEST() is required and it
 *    initializes thread-MPI with the specified number of threads and runs the
 *    rest of the test on each of the threads.
 *
 * You need to be extra careful for variables in the test fixture, if you use
 * one: when run under thread-MPI, these will be shared across all the ranks,
 * while under real MPI, these are naturally different for each process.
 * Local variables in the test body are private to each rank in both cases.
 *
 * Currently, it is not possible to specify the number of ranks as one, because
 * that will lead to problems with (at least) thread-MPI, but such tests can be
 * written as serial tests anyways.
 *
 * \ingroup module_testutils
 */
#if GMX_THREAD_MPI
#    define GMX_MPI_TEST(expectedRankCount)                                                 \
        do                                                                                  \
        {                                                                                   \
            ASSERT_EQ(expectedRankCount, ::gmx::test::getNumberOfTestMpiRanks());           \
            using MyTestClass = std::remove_reference_t<decltype(*this)>;                   \
            if (!::gmx::test::threadMpiTestRunner(std::bind(&MyTestClass::TestBody, this))) \
            {                                                                               \
                return;                                                                     \
            }                                                                               \
        } while (0)
#else
#    define GMX_MPI_TEST(expectedRankCount) \
        ASSERT_EQ(expectedRankCount, ::gmx::test::getNumberOfTestMpiRanks())
#endif

} // namespace test
} // namespace gmx

#endif
