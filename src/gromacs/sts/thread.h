/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

#ifndef STS_THREAD_H
#define STS_THREAD_H

#include <cassert>

#include <deque>
#include <memory>

#include <atomic>
#include <thread>

#include "gromacs/sts/range.h"
#include "gromacs/sts/task.h"

#if (__GNUC__ == 4 && __GNUC_MINOR__ <= 7) || (defined __ICC && __ICC <= 1400)
#define thread_local __thread
#endif
#if defined _MSC_VER && _MSC_VER == 1800
#define thread_local __declspec(thread)
#endif

/*! \brief
 * Every thread has one associated object of this type
 *
 * Both the OS created (master-) thread and all threads created by STS
 * have an object of this type. The later ones contain a std::thread
 * in thread_. Each object contains the queue of subtasks executed
 * by the respective thread during one step.
 */
class Thread {
public:
    /*! \brief
     * Constructor
     *
     * \param[in] id   Id given to this thread. 0 is the OS thread.
     */
    Thread(int id) :nextSubtaskId_(0) {
        if (id!=0) {
            thread_.reset(new std::thread([=](){id_=id; doWork();}));
        }
    }
    //! Disable move constructor
    Thread(Thread&&) = delete; //nextSubtaskId_ access would not be thread-safe
    //! Disable move assignment
    Thread& operator=(Thread&&) = delete;
    /*! \brief
     * Execute the whole queue of subtasks
     *
     * Gets executed for the OS thread by STS::wait and for STS created
     * threads from Thread::doWork
     */
    void processQueue();
    //! Execute the next subtask in the queue
    void processTask();
    //! Wait for thread to finish
    void join() { if (thread_) thread_->join(); }
    /*! \brief
     * Return thread Id
     *
     * Note: This is a static method. Id returned depends on the thread executing
     * not the object. Only use Thread::getID() not t.getId() to avoid confusion.
     *
     * \returns thread Id
     */
    static int getId() { return id_; }
    //! \brief Return current subtask Id
    int getCurrentSubTaskId() {
        return nextSubtaskId_ - 1;
    }
private:
    void doWork(); //function executed by worker threads

    unsigned int nextSubtaskId_;
    std::unique_ptr<std::thread> thread_;
    static thread_local int id_;
};

#endif // STS_THREAD_H
