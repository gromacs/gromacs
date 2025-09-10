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
/*! \internal \file
 *
 * \brief This declares a helper routine for coordinating behavior
 * when only a subset of MPI ranks might throw.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_COORDINATE_EXCEPTION_HANDLING_H
#define GMX_UTILITY_COORDINATE_EXCEPTION_HANDLING_H

#include <exception>
#include <type_traits>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{
namespace detail
{

/*! \brief If any rank threw, make sure all throw
 *
 * Any rank that threw an exception will re-throw it. If so, any other
 * rank will throw \c ParallelConsistencyError. Otherwise, does not
 * throw.
 *
 * Does global MPI communication on \c communicator.
 *
 * \param[in]  communicator   The MPI communicator over which exception-throwing
 *                            behaviour should be coordinated. Can be MPI_COMM_NULL
 *                            only for thread-MPI before rank launch.
 * \param[in]  thisRankThrew  Describes whether the current
 *                            MPI rank threw an exception
 * \param[in]  exceptionPtr   If the current MPI rank threw an exception, points
 *                            to that exception, otherwise nullptr.
 */
static void coordinateMpiRanks(MPI_Comm communicator, const bool thisRankThrew, const std::exception_ptr exceptionPtr)
{
    // There is nothing to do with no MPI, or before thread-MPI has
    // launched ranks, as there is only one rank. We default this
    // assignment to the value for the single-rank case, for
    // simplicity.
    int numRanksThatThrew = static_cast<int>(thisRankThrew);
    if (GMX_LIB_MPI or (GMX_THREAD_MPI and communicator != MPI_COMM_NULL))
    {
#if GMX_THREAD_MPI
        // Even after launching ranks, thread-MPI can't do collectives
        // on single-rank communicators.
        int numRanks;
        MPI_Comm_size(communicator, &numRanks);
        if (numRanks > 1)
#endif
        {
#if GMX_MPI
            int thisRankThrewInt = static_cast<int>(thisRankThrew);
            MPI_Allreduce(&thisRankThrewInt, &numRanksThatThrew, 1, MPI_INT, MPI_SUM, communicator);
#endif
        }
    }

    // Throw in a globally coordinated way, if needed
    if (numRanksThatThrew > 0)
    {
        if (exceptionPtr)
        {
            std::rethrow_exception(exceptionPtr);
        }
        else
        {
            GMX_THROW(ParallelConsistencyError("Another MPI rank encountered an exception"));
        }
    }
}

} // namespace detail

/*! \brief Call the \c callable and then make sure all MPI ranks of \c
 * communicator have consistent exception-throwing behaviour.
 *
 * Any rank that threw an exception in \c callable will re-throw
 * it. If so, any other rank will throw \c
 * ParallelConsistencyError. Otherwise, does not throw.
 *
 * Note that the \c callable must not require any arguments to be
 * passed, so judicious use of a function object, lambda or
 * std::function may be required.
 *
 * \param[in]  communicator  The MPI communicator over which exception-throwing
 *                           behaviour should be coordinated
 * \param[in]  callable      The "callable," which can be a function, function
 *                           object, lambda, or std::function object. It may
 *                           return a value, but is not required to do so.
 */
template<typename Callable>
auto coordinateExceptionHandling(MPI_Comm communicator, Callable callable) ->
        typename std::invoke_result_t<Callable>
{
    bool               thisRankThrew = false;
    std::exception_ptr exceptionPtr;

    if constexpr (std::is_same_v<typename std::invoke_result<Callable>::type, void>)
    {
        // Handle the simple case where there was no return value
        try
        {
            callable();
        }
        catch (const std::exception& /*ex*/)
        {
            exceptionPtr  = std::current_exception();
            thisRankThrew = true;
        }
        detail::coordinateMpiRanks(communicator, thisRankThrew, exceptionPtr);
    }
    else
    {
        // Handle the complex case, where the return value of the
        // callable must be captured before coordinating on a possible
        // thrown exception on another rank, before returning the
        // returned value when in fact no exception was thrown
        // anywhere.
        typename std::invoke_result<Callable>::type returnValues;
        try
        {
            returnValues = callable();
        }
        catch (const std::exception& /*ex*/)
        {
            exceptionPtr  = std::current_exception();
            thisRankThrew = true;
        }
        detail::coordinateMpiRanks(communicator, thisRankThrew, exceptionPtr);
        // returnValues is uninitialized only on a rank where the
        // callable threw, and in that case coordinateMpiRanks ensures
        // all ranks throw, so an uninitialized value is never
        // returned, so should not be warned about.
#if defined(__GNUC__)
        GCC_DIAGNOSTIC_IGNORE("-Wmaybe-uninitialized")
#endif
        return returnValues;
#if defined(__GNUC__)
        GCC_DIAGNOSTIC_RESET
#endif
    }
}

} // namespace gmx

#endif
