/*
   This source code file is part of thread_mpi.
   Written by Sander Pronk, Erik Lindahl, and possibly others.

   Copyright (c) 2009, Sander Pronk, Erik Lindahl.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1) Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2) Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   3) Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   If you want to redistribute modifications, please consider that
   scientific software is very special. Version control is crucial -
   bugs must be traceable. We will be happy to consider code for
   inclusion in the official distribution, but derived work should not
   be called official thread_mpi. Details are found in the README & COPYING
   files.
 */

/** \file
 *
 * \brief mutex objects with C++11 API compatibility.
 *
 * This header contains classes mutex and lock_guard, used in C++ mutex
 * implementations safe for exceptions.
 */

#ifndef TMPI_MUTEX_H_
#define TMPI_MUTEX_H_

#include "visibility.h"
#include "system_error.h"
#include "threads.h"

#ifdef __cplusplus


namespace tMPI
{
/*! \brief A lock guard class that allows for the simple management of
           mutexes. C++11 compatible.

    In C++, mutexes would normally have to be unlocked with explicit
    exception handlers and unlock statements. This class automates that
    by handling the mutex unlock in a destructor. The constructor locks
    the mutex.

    Usage example:
    tMPI::mutex mtx;
    void do_count()
    {
        tMPI::lock_guard<tMPI::mutex> lock(mtx);
        count += 1;
    }
 */
template <class Mutex> class TMPI_EXPORT lock_guard
{
    public:
        //! Lockable type that this lock operates on.
        typedef Mutex mutex_type;
        /*! \brief The constructor, which locks the mutex.

            \param m The exisiting (globally accessible) mutex to lock. */
        explicit lock_guard(mutex_type &m) : m_(m)
        {
            m_.lock();
        }
        //lock_guard(mutex_type &m, adopt_lock_t t);

        /*! \brief The destructor, which unlocks the mutex */
        ~lock_guard()
        {
            m_.unlock();
        }
    private:
        // forbid copy constructor & assignment
        lock_guard(const lock_guard &l);
        lock_guard &operator=(const lock_guard &l);

        mutex_type &m_;
};

/*! \brief A basic mutex class with C++11 compatibility.  */
class TMPI_EXPORT mutex
{
    public:
        //! Type of the native mutex handle.
        typedef tMPI_Thread_mutex_t* native_handle_type;

        /*! \brief The constructor.

           Throws a tMPI::system_error exception upon failure. */
        mutex()
        {
            int ret = tMPI_Thread_mutex_init(&handle_);
            if (ret)
            {
                throw system_error(ret);
            }
        }

        /*! \brief The destructor.*/
        ~mutex()
        {
            tMPI_Thread_mutex_destroy(&handle_);
        }

        /*! \brief The lock function.

           Throws a tMPI::system_error exception upon failure. */
        void lock()
        {
            int ret = tMPI_Thread_mutex_lock(&handle_);
            if (ret)
            {
                throw system_error(ret);
            }
        }

        /*! \brief The try_lock function.

           Throws a tMPI::system_error exception upon failure.
           \return true if the lock was locked successfully, false if not*/
        bool try_lock()
        {
            if (tMPI_Thread_mutex_trylock(&handle_))
            {
                return false;
            }
            return true;
        }

        /*! \brief The unlock function.

           Throws a tMPI::system_error exception upon failure. */
        void unlock()
        {
            int ret = tMPI_Thread_mutex_unlock(&handle_);
            if (ret)
            {
                throw system_error(ret);
            }
        }

        //! Returns the native handle for this mutex.
        native_handle_type native_handle() { return &handle_; }
    private:
        // forbid copy constructor & assignment
        mutex(const mutex &m);
        mutex &operator=(const mutex &m);

        tMPI_Thread_mutex_t handle_;
};
}

#endif /* __cplusplus */

#endif /* TMPI_MUTEX_H_ */
