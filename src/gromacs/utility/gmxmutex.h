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
#ifndef GMX_UTILITY_GMXMUTEX_H
#define GMX_UTILITY_GMXMUTEX_H

/*! \file
 * \brief
 * Defines mutex class for use in \Gromacs.
 *
 * \todo Consider this with std::mutex. There are currently no per-MD
 * step mutex operations. All uses are either wrapping I/O open/close
 * operations, or as part of setup and cleanup, so performance is of
 * low importance.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */

#include "thread_mpi/threads.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! This class implements a mutex suitable for RAII-style use.
 *
 * The set of operations are modelled on std::mutex, to which we might
 * switch some time. Thus, this class meets the requirements of
 * BasicLockable and Lockable, so is appropriate to use with
 * std::unique_lock and std::lock_guard. The use of such types is
 * recommended over the use of raw Mutex objects, just as for
 * std::mutex.
 *
 * This class is implemented by thread-MPI functionality.
 */
class Mutex
{
    public:
        //! Constructor.
        Mutex() : mutex_(TMPI_THREAD_MUTEX_INITIALIZER)
        {
            tMPI_Thread_mutex_init(&mutex_);
        }
        Mutex(const Mutex &)    = delete;
        Mutex &operator=(Mutex) = delete;
        //! Destructor.
        ~Mutex()
        {
            tMPI_Thread_mutex_destroy(&mutex_);
        }
        /*! \brief Locks the mutex, once it becomes available.
         *
         * Since a mutex should generally not be used in
         * performance-critical code, the additional overhead of
         * checking and handling the error condition is of no concern.
         *
         * \throws InternalError  If a mutex cannot be obtained. This mirrors
         *     the behaviour of std::mutex.
         */
        void lock()
        {
            int ret;
            if ((ret = tMPI_Thread_mutex_lock(&mutex_)) != 0)
            {
                GMX_THROW(InternalError(formatString("Failed to lock mutex. Return value was %d.", ret)));
            }
        }
        //! Try to obtain the mutex if it is available, else return false.
        bool try_lock()
        {
            return tMPI_Thread_mutex_trylock(&mutex_) == 0;
        }
        //! Release the mutex.
        void unlock()
        {
            tMPI_Thread_mutex_unlock(&mutex_);
        }
    private:
        tMPI_Thread_mutex_t mutex_;
};

} // namespace

#endif
