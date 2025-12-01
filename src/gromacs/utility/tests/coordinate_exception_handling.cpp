/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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

#include "gromacs/utility/coordinate_exception_handling.h"

#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include "testutils/mpitest.h"

namespace gmx
{
namespace test
{
namespace
{

void functionReturningVoid() {}

struct FunctionObjectReturningVoid
{
    void operator()() const {}
};

TEST(CoordinateExceptionHandlingTest, CallablesWithoutReturnValueNeverThrow)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    EXPECT_NO_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, functionReturningVoid));
    EXPECT_NO_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, FunctionObjectReturningVoid{}));
    const auto testLambda = []() { return; };
    EXPECT_NO_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testLambda));
    const std::function testStdFunction{ []() { return; } };
    EXPECT_NO_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testStdFunction));
}

} // namespace

namespace
{

int functionReturningInt()
{
    return 1;
}

std::unique_ptr<int> functionReturningUniquePtr()
{
    return std::make_unique<int>(1);
}

struct FunctionObjectReturningInt
{
    int operator()() const { return 1; }
};

TEST(CoordinateExceptionHandlingTest, CallablesWithReturnValueAlwaysReturnThatValue)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    EXPECT_EQ(coordinateExceptionHandling(MPI_COMM_WORLD, functionReturningInt), 1);
    EXPECT_EQ(*coordinateExceptionHandling(MPI_COMM_WORLD, functionReturningUniquePtr), 1)
            << "returning move-only types works";
    EXPECT_EQ(coordinateExceptionHandling(MPI_COMM_WORLD, FunctionObjectReturningInt{}), 1);
    const auto testLambda = []() -> int { return 1; };
    EXPECT_EQ(coordinateExceptionHandling(MPI_COMM_WORLD, testLambda), 1);
    const std::function testStdFunction{ []() -> int { return 1; } };
    EXPECT_EQ(coordinateExceptionHandling(MPI_COMM_WORLD, testStdFunction), 1);
}

} // namespace

namespace
{

[[noreturn]] void functionReturningVoidThatThrows()
{
    throw std::exception();
}

struct FunctionObjectReturningVoidThatThrows
{
    [[noreturn]] void operator()() const { throw std::exception(); }
};

TEST(CoordinateExceptionHandlingTest, CallablesWithoutReturnValueThatThrowAlwaysThrow)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, functionReturningVoidThatThrows),
                 std::exception);
    EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, FunctionObjectReturningVoidThatThrows{}),
                 std::exception);
    // oneAPI 2025.3 warns about the absence of a noreturn attribute
    // when compiling for C++17, but actually using an attribute is
    // not supported until C++23. Fortunately we just don't care about
    // this for a test case.
    CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn");
    const auto testLambda = []() { throw std::exception(); };
    EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testLambda), std::exception);
    const std::function f{ []() { throw std::exception(); } };
    EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, f), std::exception);
    CLANG_DIAGNOSTIC_RESET;
}

} // namespace

namespace
{

void functionReturningVoidThatThrowsOnOneRank()
{
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if (mpiRank == 0)
    {
        throw std::exception();
    }
}

struct FunctionObjectReturningVoidThatThrowsOnOneRank
{
    FunctionObjectReturningVoidThatThrowsOnOneRank(const int mpiRank) : mpiRank_(mpiRank) {}
    void operator()() const
    {
        if (mpiRank_ == 0)
        {
            throw std::exception();
        }
    }

private:
    const int mpiRank_;
};

TEST(CoordinateExceptionHandlingTest, CallablesWithoutReturnValueThatThrowsOnOneRankAlwaysThrowCorrectly)
{
    GMX_MPI_TEST(RequireMinimumRankCount<2>);
    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    const auto testLambda = [mpiRank]()
    {
        if (mpiRank == 0)
        {
            throw std::exception();
        }
    };
    const std::function testStdFunction{ [mpiRank]()
                                         {
                                             if (mpiRank == 0)
                                             {
                                                 throw std::exception();
                                             }
                                         } };
    if (mpiRank == 0)
    {
        EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, functionReturningVoidThatThrowsOnOneRank),
                     std::exception);
        EXPECT_THROW(coordinateExceptionHandling(
                             MPI_COMM_WORLD, FunctionObjectReturningVoidThatThrowsOnOneRank{ mpiRank }),
                     std::exception);
        EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testLambda), std::exception);
        EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testStdFunction), std::exception);
    }
    else
    {
        EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, functionReturningVoidThatThrowsOnOneRank),
                     gmx::ParallelConsistencyError);
        EXPECT_THROW(coordinateExceptionHandling(
                             MPI_COMM_WORLD, FunctionObjectReturningVoidThatThrowsOnOneRank{ mpiRank }),
                     gmx::ParallelConsistencyError);
        EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testLambda), gmx::ParallelConsistencyError);
        EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testStdFunction),
                     gmx::ParallelConsistencyError);
    }
}

} // namespace

namespace
{

TEST(CoordinateExceptionHandlingTest, NulledCallablesAlwaysThrow)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    const std::function<void(void)> testStdFunction = nullptr;
    EXPECT_THROW(coordinateExceptionHandling(MPI_COMM_WORLD, testStdFunction), std::bad_function_call);
}


} // namespace
} // namespace test
} // namespace gmx
