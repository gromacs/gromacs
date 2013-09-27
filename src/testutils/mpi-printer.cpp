/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include "mpi-printer.h"

#ifdef GMX_LIB_MPI

#include "config.h"
#include "gromacs/utility/gmxmpi.h"

#include <boost/scoped_ptr.hpp>
#include <gtest/gtest.h>

#include "gromacs/utility/common.h"
#include "gromacs/utility/gmxassert.h"

#define FORWARD_TO_DEFAULT_PRINTER1(MethodName, Param1Type) \
    virtual void MethodName(const Param1Type &param1) { \
        if (rank_ == 0) { \
            defaultPrinter_->MethodName(param1); \
        } \
    }
#define FORWARD_TO_DEFAULT_PRINTER2(MethodName, Param1Type, Param2Type) \
    virtual void MethodName(const Param1Type &param1, Param2Type param2) { \
        if (rank_ == 0) { \
            defaultPrinter_->MethodName(param1, param2); \
        } \
    }

namespace
{

class MPIEventForward : public ::testing::TestEventListener
{
    public:
        MPIEventForward(TestEventListener* defaultPrinter, int rank, int size)
            : defaultPrinter_(defaultPrinter), rank_(rank), size_(size)
        {
        }

        FORWARD_TO_DEFAULT_PRINTER1(OnTestProgramStart, ::testing::UnitTest);
        FORWARD_TO_DEFAULT_PRINTER2(OnTestIterationStart, ::testing::UnitTest, int);
        FORWARD_TO_DEFAULT_PRINTER1(OnEnvironmentsSetUpStart, ::testing::UnitTest);
        FORWARD_TO_DEFAULT_PRINTER1(OnEnvironmentsSetUpEnd, ::testing::UnitTest);
        FORWARD_TO_DEFAULT_PRINTER1(OnTestCaseStart, ::testing::TestCase);
        FORWARD_TO_DEFAULT_PRINTER1(OnTestStart, ::testing::TestInfo);
        virtual void OnTestPartResult(const ::testing::TestPartResult & /*result*/)
        {
            // Do nothing; all printing is done in OnTestEnd().
        }
        virtual void OnTestEnd(const ::testing::TestInfo &test_info);
        FORWARD_TO_DEFAULT_PRINTER1(OnTestCaseEnd, ::testing::TestCase);
        FORWARD_TO_DEFAULT_PRINTER1(OnEnvironmentsTearDownStart, ::testing::UnitTest);
        FORWARD_TO_DEFAULT_PRINTER1(OnEnvironmentsTearDownEnd, ::testing::UnitTest);
        FORWARD_TO_DEFAULT_PRINTER2(OnTestIterationEnd, ::testing::UnitTest, int);
        FORWARD_TO_DEFAULT_PRINTER1(OnTestProgramEnd, ::testing::UnitTest);

    private:
        boost::scoped_ptr<TestEventListener> defaultPrinter_;
        int                                  rank_;
        int                                  size_;

        GMX_DISALLOW_COPY_AND_ASSIGN(MPIEventForward);
};

void MPIEventForward::OnTestEnd(const ::testing::TestInfo &test_info)
{
    // Serialize printing test results to stdout in rank order by
    // passing a flag in order from ranks 0 .. (n-1). Rank 0 does not
    // need to recieve, rank n-1 does not need to send.
    int timeToPrint = (0 == rank_);
    if (!timeToPrint)
    {
        MPI_Recv(&timeToPrint, 1, MPI_INT, rank_ - 1, 0, MPI_STATUS_IGNORE);
    }

    // Now this rank can print
    int local_passed = true;
    const ::testing::TestResult *result = test_info.result();
    if (result->Failed())
    {
        printf("Test failures from rank %d:\n", rank_);
        for (int i = 0; i < result->total_part_count(); ++i)
        {
            defaultPrinter_->OnTestPartResult(result->GetTestPartResult(i));
        }
        local_passed = false;
    }

    // Pass on the printing token to the next rank
    if (size_ != rank_ + 1)
    {
        MPI_Send(&timeToPrint, 1, MPI_INT, rank_ + 1, 0, MPI_STATUS_IGNORE);
    }

    int all_passed = true;
    MPI_Reduce(&local_passed, &all_passed, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);

    if (rank_ == 0)
    {
        if (!all_passed && local_passed)
        {
            // This marks the current test failed, and modifies test_info
            // behind the scenes, so the default printer sees the test as
            // failed even if local_passed is true.
            ADD_FAILURE();
        }

        defaultPrinter_->OnTestEnd(test_info);
    }
}

} // namespace

#endif

/*! \brief Deploy custom GTest event listeners for handling errors
 * when run under MPI.
 *
 * Only one rank should report the test result. Errors detected on a
 * subset of ranks need to be reported individually, and as an overall
 * failure. */
void gmx::test::initMPIOutput()
{
#ifdef GMX_LIB_MPI
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        return;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ::testing::UnitTest           &unit_test  = *::testing::UnitTest::GetInstance();
    ::testing::TestEventListeners &listeners  = unit_test.listeners();
    ::testing::TestEventListener  *defprinter =
        listeners.Release(listeners.default_result_printer());
    listeners.Append(new MPIEventForward(defprinter, rank, size));
#endif
}
