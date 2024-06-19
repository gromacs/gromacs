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
/*! \internal \file
 * \brief
 * Implements functions in mpitest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/mpitest.h"

#include "config.h"

#include <functional>

#include <gtest/gtest.h>

#include "thread_mpi/tmpi.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{

#if GMX_THREAD_MPI

namespace
{

//! Number of tMPI threads.
int g_numThreads = 1;
//! \cond
GMX_TEST_OPTIONS(ThreadMpiTestOptions, options)
{
    options->addOption(
            IntegerOption("ntmpi").store(&g_numThreads).description("Number of thread-MPI threads/ranks for the test"));
}
//! \endcond

//! Thread entry function for other thread-MPI threads.
void threadStartFunc(const void* data)
{
    const std::function<void()>& testBody = *reinterpret_cast<const std::function<void()>*>(data);
    try
    {
        testBody();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

//! Helper function for starting thread-MPI threads for a test.
bool startThreads(std::function<void()>* testBody)
{
    int ret = tMPI_Init_fn(TRUE, g_numThreads, TMPI_AFFINITY_NONE, threadStartFunc, testBody);
    return ret == TMPI_SUCCESS;
}

class InTestGuard
{
public:
    explicit InTestGuard(bool* inTest) : inTest_(inTest) { *inTest = true; }
    ~InTestGuard() { *inTest_ = false; }

private:
    bool* inTest_;
};

} // namespace

//! \cond internal
bool threadMpiTestRunner(std::function<void()> testBody)
{
    static bool inTest = false;

    if (inTest || g_numThreads <= 1)
    {
        return true;
    }
#    if GMX_THREAD_MPI && !defined(GTEST_IS_THREADSAFE)
    ADD_FAILURE() << "Google Test is not thread safe on this platform. "
                  << "Cannot run multi-rank tests with thread-MPI.";
#    else
    InTestGuard guard(&inTest);
    if (!startThreads(&testBody))
    {
        return false;
    }
    try
    {
        testBody();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    tMPI_Finalize();
#    endif
    return false;
}
//! \endcond

#endif

int getNumberOfTestMpiRanks()
{
#if GMX_THREAD_MPI
    return g_numThreads;
#else
    return gmx_node_num();
#endif
}

} // namespace test
} // namespace gmx
