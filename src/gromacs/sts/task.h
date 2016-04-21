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

#ifndef STS_TASK_H
#define STS_TASK_H

#include <map>
#include <vector>

#include <chrono>

#include "gromacs/sts/barrier.h"
#include "gromacs/sts/range.h"

using sts_clock = std::chrono::steady_clock;

//! \internal Interface of the executable function of a task
class ITaskFunctor {
public:
    /*! \brief
     * Run the function of the task
     *
     * \param[in] range  range of task to be executed. Ignored for basic task.
     */
    virtual void run(Range<Ratio> range) = 0;
    virtual ~ITaskFunctor() {};
};

//! \internal Loop task functor
template<class F>
class LoopTaskFunctor : public ITaskFunctor {
public:
    /*! \brief
     * Constructor
     *
     * \param[in] f    lambda of loop body
     * \param[in] r    range of loop
     */
    LoopTaskFunctor<F>(F f, Range<int64_t> r): body_(f), range_(r) {}
    void run(Range<Ratio> r) {
        Range<int64_t> s = range_.subset(r); //compute sub-range of this execution
        for (int i=s.start; i<s.end; i++) {
            body_(i);
        }
    }
private:
    F body_;
    Range<int64_t> range_;
};

//! \internal Basic (non-loop) task functor
template<class F>
class BasicTaskFunctor : public ITaskFunctor {
public:
    /*! \brief
     * Constructor
     *
     * \param[in] f    lambda of function
     */
    BasicTaskFunctor<F>(F f) : func_(f) {};
    void run(Range<Ratio>) {
        func_();
    }
private:
    F func_;
};

/*! \internal \brief
 * The part of a task done by one thread
 *
 * Contains the input to a thread what to do (taskId_, range_) and the output
 * from the thread upon completion of the sub-task (done_, timing).
 * For a loop task the range_ is the subsection done by a thread. A basic
 * task is always completely executed by a single thread.
 */
class SubTask {
public:
    /*! \brief
     * Constructor
     *
     * \param[in] taskId   The ID of the task this is part of.
     * \param[in] range    Out of a possible range from 0 to 1, the section in
                           this part. Ignored for basic tasks.
     */
    SubTask(int taskId, Range<Ratio> range) : taskId_(taskId), range_(range) {}
    int getTaskId() const { return taskId_; }               //!< get Task Id
    const Range<Ratio>& getRange() const { return range_; } //!< get Range

    sts_clock::duration waitTime_; /**< Time spent until task was ready */
    sts_clock::duration runTime_;  /**< Time spent executing sub-task  */
private:
    int taskId_;                   /**< The ID of the task this is a part of */
    Range<Ratio> range_;           /**< Range (out of [0,1]) of loop part */
};

/*! \internal \brief
 * A task to be executed
 *
 * Can either be a function or loop. Depending on the schedule is
 * executed synchronous or asynchronous. Functions are always
 * executed by a single thread. Loops are executed, depending on
 * the schedule, in serial or in parallel.
 */
struct Task {
    void *reduction_; //!< Reduction function to execute after task completes
    ITaskFunctor *functor_;      //!< The function/loop to execute
    MOBarrier functorBeginBarrier_; //!< Many-to-one barrier to sync threads at beginning of loop
    OMBarrier functorEndBarrier_; //!< One-to-many barrier to sync threads at end of loop
    //! All subtasks of this task. One for each section of a loop. One for a basic task.
    std::vector<SubTask const*> subtasks_; //!< Subtasks to be executed by a single thread
    //!< The waiting time in the implied barrier at the end of a loop. Zero for basic task.
    sts_clock::duration waitTime_; //!< Time that main thread waits on end barrier
    sts_clock::duration reductionTime_; //!< Time spent doing reduction

    Task() :reduction_(nullptr), functor_(nullptr), waitTime_(0), reductionTime_(0), numThreads_(0) {}
    /*! \brief
     * Add a new subtask for this task
     *
     * \param[in] threadId  thread to which this subtask is assigned
     * \param[in] t         subtask
     */
    void pushSubtask(int threadId, SubTask const* t) {
        subtasks_.push_back(t);
        if (threadTaskIds_.find(threadId) == threadTaskIds_.end()) {
            threadTaskIds_[threadId] = numThreads_;
            numThreads_++;
        }
    }
    //! \brief Remove all subtasks for this task
    void clearSubtasks() {
        subtasks_.clear();
        threadTaskIds_.clear();
        numThreads_ = 0;
    }
    //! \brief Get total number of threads assigned to some subtask for this task
    int getNumThreads() const {
        return numThreads_;
    }
    /*! \brief
     * Get task-specific thread Id for the given STS thread Id
     *
     * \param[in] threadId  STS thread id
     */
    int getThreadId(int threadId) const {
        auto id = threadTaskIds_.find(threadId);
        if (id == threadTaskIds_.end()) {
            return -1;
        }
        return (*id).first;
    }
private:
    int numThreads_;
    //! Map STS thread id to an id only for this task (task ids are consecutive starting from 0)
    std::map<int, int> threadTaskIds_;
};

#endif // STS_TASK_H
