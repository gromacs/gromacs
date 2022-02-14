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
#include <string>
#include <type_traits>

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

/*! \brief Implementation of MPI test runner for thread-MPI
 *
 * See documentation GMX_MPI_TEST */
#if GMX_THREAD_MPI
#    define GMX_MPI_TEST_INNER                                                                  \
        do                                                                                      \
        {                                                                                       \
            using MyTestClass = std::remove_reference_t<decltype(*this)>;                       \
            if (!::gmx::test::threadMpiTestRunner([this]() { this->MyTestClass::TestBody(); })) \
            {                                                                                   \
                return;                                                                         \
            }                                                                                   \
        } while (0)
#else
#    define GMX_MPI_TEST_INNER
#endif

/*! \brief Declares that this test is an MPI-enabled unit test and
 * expresses the conditions under which it can run.
 *
 * To write unit tests that run under MPI, you need to do a few things:
 *  - Put GMX_MPI_TEST(RankRequirement) as the first statement in your
 *    test body and either declare or use a suitable class as its
 *    argument to express what requirements exist on the number of MPI
 *    ranks for this test.
 *  - Declare your unit test in CMake with gmx_add_mpi_unit_test().
 *    Note that all tests in the binary should fulfill the conditions above.
 *
 * When you do the above, the following will happen:
 *  - The test will get compiled only if thread-MPI or real MPI is enabled.
 *  - The test will get executed only when the specified condition on
 *    the the number of ranks is satisfied.
 *  - If you are using real MPI, the whole test binary is run under
 *    MPI and test execution across the processes is synchronized
 *    (GMX_MPI_TEST() actually has no effect in this case, the
 *    synchronization is handled at a higher level).
 *  - If you are using thread-MPI, GMX_MPI_TEST() is required and it
 *    initializes thread-MPI with the specified number of threads and
 *    runs the rest of the test on each of the threads.
 *
 * \param[in] RankRequirement Class that expresses the necessary
 *     conditions on the number of MPI ranks for the test to continue.
 *     If run with unsupported number of ranks, the remainder of the
 *     test body is skipped, and the GTEST_SKIP() mechanism used to
 *     report the reason why the number of MPI ranks is unsuitable.
 *
 * The RankRequirement class must have two static members; a static
 * method \c bool conditionSatisfied(const int) that can be passed the
 * number of ranks present at run time and return whether the test can
 * run with that number of ranks, and a static const string \c
 * s_skipReason describing the reason why the test cannot be run, when
 * that is the case.
 *
 * You need to be extra careful for variables in the test fixture, if you use
 * one: when run under thread-MPI, these will be shared across all the ranks,
 * while under real MPI, these are naturally different for each process.
 * Local variables in the test body are private to each rank in both cases.
 *
 * Currently, it is not possible to require the use of a single MPI
 * rank, because that will lead to problems with (at least)
 * thread-MPI, but such tests can be written as serial tests anyway.
 *
 * \ingroup module_testutils
 */
#define GMX_MPI_TEST(RankRequirement)                                                         \
    const int numRanks = ::gmx::test::getNumberOfTestMpiRanks();                              \
    if (!RankRequirement::conditionSatisfied(numRanks))                                       \
    {                                                                                         \
        GTEST_SKIP() << std::string("Test skipped because ") + RankRequirement::s_skipReason; \
        return;                                                                               \
    }                                                                                         \
    GMX_MPI_TEST_INNER;

//! \internal \brief Helper for GMX_MPI_TEST to permit any rank count
class AllowAnyRankCount
{
public:
    /*! \brief Function called by GMX_MPI_CONDITIONAL_TEST to see
     * whether the test conditions are satisifed */
    static bool conditionSatisfied(const int /* numRanks */) { return true; }
    //! Reason to echo when skipping the test
    inline static const char* s_skipReason = "UNUSED - any rank count satisfies";
};

//! \internal \brief Helper for GMX_MPI_TEST to permit only a specific rank count
template<int requiredNumRanks>
class RequireRankCount
{
public:
    //! Function to require a specific number of ranks
    static bool conditionSatisfied(const int numRanks) { return numRanks == requiredNumRanks; }
    //! Text to echo when skipping a test that does not satisfy the requirement
    inline static const std::string s_skipReason =
            std::to_string(requiredNumRanks) + " ranks are required";
};

//! \internal \brief Helper for GMX_MPI_TEST to permit only a specific rank count
template<int minimumNumRanks>
class RequireMinimumRankCount
{
public:
    //! Function to require at least the minimum number of ranks
    static bool conditionSatisfied(const int numRanks) { return numRanks >= minimumNumRanks; }
    //! Text to echo when skipping a test that does not satisfy the requirement
    inline static const std::string s_skipReason =
            std::to_string(minimumNumRanks) + " or more ranks are required";
};


} // namespace test
} // namespace gmx

#endif
