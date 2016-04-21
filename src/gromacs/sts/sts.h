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

#ifndef STS_STS_H
#define STS_STS_H

#include <cassert>

#include <deque>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <atomic>
#include <chrono>

#include "gromacs/sts/range.h"
#include "gromacs/sts/reduce.h"
#include "gromacs/sts/task.h"
#include "gromacs/sts/thread.h"

/* Overall design:
 * The framework can execute simple tasks (via "run") and execute loops in
 * parallel (via "parallel_for"). It supports two run modi: either with an
 * explicit schedule or with a default schedule. With the default schedule
 * tasks are run in serial and only loop level parallelism is used. This is
 * useful if either the tasks are not yet known or only simple parallelism is
 * needed. With an explicit schedule one can specify which task runs on which
 * thread and in which order (based on the order of task assignment). Loops
 * can be split among the threads using ratios (e.g. thread 0 does 1/3 of
 * the loop while thread 1 does the remaining 2/3). The idea is that this
 * schedule is either computed by the user of the framework using "assign"
 * or automatically computed by the framework using "reschedule." (Automatic
 * scheduling is not yet implemented.) Timing data is recorded for each task
 * so that adjustments can be made (or not) after each "step." One "step"
 * contains a number of scheduled tasks and a new step starts when "nextStep"
 * is called. Normally, a step will be one iteration of a main loop, like a
 * time step in MD, but this is of course not required. The part of a task
 * done by a thread is called a sub-task. A simple task is always fully
 * done by one thread and for a loop-task the range done by each thread is
 * specified. The whole design is lock free and only relies on atomics.
 */

/*! \brief
 * Static task scheduler
 *
 * Allows running an asynchronous function with run() and execute loops in parallel
 * with parallel_for(). The default schedule only uses loop level parallelism and
 * executes run() functions synchronously. A schedule with task level parallelism
 * is created automatically by calling reschedule() or manual by calling assign().
 * After the queue of tasks in one schedule are completed, one can do one of
 * three things for the next step:
 * a) call nextStep() and reuse the schedule
 * b) call reschedule() to let the scheduler automatically compute a new schedule
 * c) call clearAssignments(), assign(), and nextStep() to manual specify schedule
 */
class STS {
public:
    const std::string id;
    /*! \brief
     * Startup STS and set the number of threads.
     * No STS functions should be called before Startup.
     */
    static void startup(int numThreads) {
        assert(numThreads > 0);
        // Allow multiple calls but make sure thread count is the same for all.
        if (threads_.size() > 0) {
            assert(numThreads == threads_.size());
            return;
        }
        for (int id = 0; id < numThreads; id++) {
            threads_.emplace_back(id); //create threads
        }
        // Create the default STS instance, which uses a default schedule.
        // This schedule provides a quick way to parallelize loops.
        defaultInstance_ = new STS("default");
        defaultInstance_->setDefaultSchedule();
        instance_ = defaultInstance_;
    }
    /*! \brief
     * Stops all threads.
     * No STS functions should be called after Shutdown.
     */
    static void shutdown() {
        assert(instance_ == defaultInstance_);
        //-1 notifies threads to finish
        stepCounter_.store(-1, std::memory_order_release);
        for (unsigned int i=1;i<threads_.size();i++) {
            threads_[i].join();
        }
        delete defaultInstance_;
    }
    STS(std::string name = "") :id(name) {
        int n = getNumThreads();
        assert(n > 0);
        threadSubTasks_.resize(n);
        if (!id.empty()) {
            stsInstances_[id] = this;
        }
    }
    ~STS() {
        for (auto& taskList : threadSubTasks_) {
            for (SubTask *t : taskList) {
                delete t;
            }
        }
        if (!id.empty()) {
            stsInstances_.erase(id);
        }
    }
    /*! \brief
     * Get number of threads in the pool
     */
    static int getNumThreads() {
        auto n = threads_.size();
        assert(n > 0);
        return n;
    }
    /*! \brief
     * Assign task to a thread
     *
     * If a range for a loop task is specified, only that section of the loop is assigned.
     * In that case it is important to assign the remaining loop out of [0,1] also to
     * some other thread. It is valid to assign multiple parts of a loop to the same thread.
     * The order of assign calls specifies in which order the thread executes the tasks.
     *
     * \param[in] label    The label of the task. Needs to match the run()/parallel_for() label
     * \param[in] threadId The Id of the thread to assign to
     * \param[in] range    The range for a loop task to assing. Ignored for basic task.
     */
    void assign(std::string label, int threadId, Range<Ratio> range = Range<Ratio>(1)) {
        int id = getTaskId(label);
        assert(range.start>=0 && range.end<=1);
        SubTask *t = new SubTask(id, range);
        threadSubTasks_[threadId].push_back(t);
        tasks_[id].pushSubtask(threadId, t);
    }
    //! Clear all assignments
    void clearAssignments() {
        for (auto &taskList : threadSubTasks_) {
            taskList.clear();
        }
        for (auto &task : tasks_) {
            task.clearSubtasks();
        }
    }
    void setDefaultSchedule() {
        bUseDefaultSchedule_ = true;
        clearAssignments();
        int numThreads = getNumThreads();
        for (int id = 0; id < numThreads; id++) {
            assign("default", id, {{id, numThreads}, {id + 1, numThreads}});
        }
    }
    /*! \brief
     * Run an asynchronous function
     *
     * \param[in] label     The task label (needs to match assign())
     * \param[in] function  The function (or lambda) to execute
     */
    template<typename F>
    void run(std::string label, F function) {
        if (!bUseDefaultSchedule_) {
            assert(this == instance_);
            // Cannot invoke an instance that is inactive
            // (nextStep has not been called)
            assert(instance_->isActive_ == true);
        }
        else {
            // Instances with a default schedule can be run at any time except
            // in the middle of an active schedule.
            assert(instance_->isActive_ == false);
        }
        if (!isTaskAssigned(label) || bUseDefaultSchedule_) {
            function();
        } else {
            tasks_[getTaskId(label)].functor_ = new BasicTaskFunctor<F>(function);
            tasks_[getTaskId(label)].functorBeginBarrier_.open();
        }
    }
    //! Notify threads to start computing the next step
    void nextStep() {
        if (!bUseDefaultSchedule_) {
            nextStepInternal();
        }
    }
    /*! \brief
     * Execute a parallel for loop
     *
     * \param[in] label    The task label (needs to match assign())
     * \param[in] start    The start index of the loop
     * \param[in] end      The end index of the loop
     * \param[in] body     The function (or lambda) to execute as loop body
     * \param[in] red      Optional reduction
     */
    template<typename F, typename T=int>
    void parallel_for(std::string label, int64_t start, int64_t end, F body, TaskReduction<T> *red = nullptr) {
        if (!bUseDefaultSchedule_) {
            assert(this == instance_);
            // Cannot invoke an instance that is inactive
            // (nextStep has not been called)
            assert(instance_->isActive_ == true);
        }
        else {
            // Instances with a default schedule can be run at any time except
            // in the middle of an active schedule.
            assert(instance_->isActive_ == false);
        }
        int taskId = -1;
        if (bUseDefaultSchedule_) {
            assert(isTaskAssigned("default"));
            taskId = getTaskId("default");
            assert(taskId == 0); // Default schedule should have only the "default" task with id 0.
            nextStepInternal();  // Default schedule has only a single step and the user doesn't need to call nextStep
        } else if (!isTaskAssigned(label)) {
            for (int i=start; i<end; i++) {
                body(i);
            }
            if (red != nullptr) {
                red->reduce();
            }
            return;
        } else {
            taskId = getTaskId(label);
        }
        auto &task = tasks_[taskId];
        task.reduction_ = red;
        task.functor_ = new LoopTaskFunctor<F>(body, {start, end});
        task.functorBeginBarrier_.open();
        auto &thread = threads_[Thread::getId()];
        //Calling processTask implies that the thread calling parallel_for participates in the loop and executes it next in queue
        assert(getSubTask(Thread::getId(), thread.getCurrentSubTaskId() + 1)->getTaskId() == taskId);
        thread.processTask();
        auto startWaitTime = sts_clock::now();
        task.functorEndBarrier_.wait();
        task.waitTime_ = sts_clock::now() - startWaitTime;
        // TODO: A smarter reduction would take place before the above wait.
        if (task.reduction_ != nullptr) {
            auto startReductionTime = sts_clock::now();
            static_cast< TaskReduction<T> *>(task.reduction_)->reduce();
            task.reductionTime_ = sts_clock::now() - startReductionTime;
        }
        // User does not need to call wait for default scheduling
        if (bUseDefaultSchedule_) {
            waitInternal();
        }
    }
    void skip_run(std::string label) {
        run(label, []{});
    }
    void skip_loop(std::string label) {
        // Return i to avoid compiler warnings about an unused parameter
        parallel_for(label, 0, 0, [](int i){return i;});
    }
    //! Automatically compute new schedule based on previous step timing
    void reschedule() {
        // not yet available
    }
    //! Wait on all tasks to finish
    void wait() {
        if (!bUseDefaultSchedule_) {
            waitInternal();
        }
    }
    /*! \brief
     * Returns STS instance for a given id or default instance if not found
     *
     * \param[in] STS instance id
     * \returns STS instance
     */
    static STS *getInstance(std::string id) {
        auto entry = stsInstances_.find(id);
        if (entry == stsInstances_.end()) {
            return defaultInstance_;
        }
        else {
            return entry->second;
        }
    }
    /*! \brief
     * Returns current STS instance
     *
     * WARNING: meant only  for internal use. Applications should use
     * "getInstance" for better error checking and clarity when using
     * multiple STS instances.
     *
     * \returns current STS instance
     */
    static STS *getCurrentInstance() {
        return instance_;
    }
    /*! \brief
     * Returns the task functor for a given task Id
     *
     * Waits on functor to be ready if the corresponding run()/parallel_for() hasn't been executed yet.
     *
     * \param[in] task Id
     * \returns task functor
     */
    ITaskFunctor *getTaskFunctor(int taskId) {
        tasks_[taskId].functorBeginBarrier_.wait();
        return tasks_[taskId].functor_;
    }
    SubTask *getSubTask(int threadId, int subTaskId) const {
        return threadSubTasks_[threadId][subTaskId];
    }
    int getNumSubTasks(int threadId) const {
        return threadSubTasks_[threadId].size();
    }
    void markSubtaskComplete(int taskId) {
        tasks_[taskId].functorEndBarrier_.markArrival();
    }
    /* \brief
     * Get number of threads for current task or 0 if no current task
     *
     * \return number of threads or 0 if no current task
     */
    int getTaskNumThreads() {
        const Task *t = getCurrentTask();
        if (t == nullptr) {
            return 0;
        }
        return t->getNumThreads();
    }
    /* \brief
     * Get number of threads for a given task
     *
     * \return number of threads
     */
    int getTaskNumThreads(std::string label) {
        assert(isTaskAssigned(label));
        int taskId = getTaskId(label);
        return tasks_[taskId].getNumThreads();
    }
    /* \brief
     * Get thread's id for its current task or -1
     *
     * \return thread task id or -1 if no current task
     */
    int getTaskThreadId() {
        const Task *t = getCurrentTask();
        if (t == nullptr) {
            return -1;
        }
        int ttid = t->getThreadId(Thread::getId());
        // Would mean that thread is currently running a task it was never assigned.
        assert(ttid > -1);
        return ttid;
    }
    /* \brief
     * Load atomic step counter
     *
     * \returns step counter
     */
    static int loadStepCounter() { return stepCounter_.load(std::memory_order_acquire); }
    /* \brief
     * Wait on atomic step counter to change
     *
     * param[in] c   last step processed by thread
     */
    static int waitOnStepCounter(int c) {return wait_until_not(stepCounter_, c);}

    /* Task reduction functions */

    /*! \brief
     * Create a TaskReduction object
     *
     * \param[in] taskName  Name of task to which reduction applies
     * \param[in] init      Initial value
     */
    template<typename T>
    TaskReduction<T> createTaskReduction(std::string taskName, T init) {
        int numThreads = getTaskNumThreads(taskName);
        return TaskReduction<T>(init, numThreads);
    }
    /*! \brief
     * Collect a value for a task's reduction. Must be called within a task.
     *
     * \param[in] a Value to be collected
     */
    template<typename T>
    void collect(T a) {
        collect(a, getTaskThreadId());
    }
private:
    // Helper functions for operations that need the current task
    int getCurrentTaskId() {
        int threadId = Thread::getId();
        int subTaskId = threads_[threadId].getCurrentSubTaskId();
        if (subTaskId == -1) {
            return -1;
        }
        return threadSubTasks_[threadId][subTaskId]->getTaskId();
    }
    const Task *getCurrentTask() {
        int taskId = getCurrentTaskId();
        if (taskId == -1) {
            return nullptr;
        }
        return &tasks_[taskId];
    }
    bool isTaskAssigned(std::string label) const {
        return (taskLabels_.find(label) != taskLabels_.end());
    }
    // Creates new ID for unknown label.
    // Creating IDs isn't thread safe. OK because assignments and run/parallel_for (if run without pre-assignment) are executed by master thread while other threads wait on nextStep.
    int getTaskId(std::string label) {
        auto it = taskLabels_.find(label);
        if (it != taskLabels_.end()) {
            return it->second;
        } else {
            assert(Thread::getId()==0); //creating thread should only be done by master thread
            unsigned int v = taskLabels_.size();
            assert(v==tasks_.size());
            tasks_.resize(v+1);
            taskLabels_[label] = v;
            return v;
        }
    }
    /* \brief
     * Collect the given value for the current task for later reduction
     *
     * \param[in] value to collect
     */
    template<typename T>
    void collect(T a, int ttid) {
        const Task *t = getCurrentTask();
        // TODO: This is a user error - calling collect outside of a task.
        // Currently, we simply ignore the call. How should it be handled?
        if (t == nullptr) {
            return;
        }
        (static_cast<TaskReduction<T> *>(t->reduction_))->collect(a, ttid);
    }
    //! Notify threads to start computing the next step
    void nextStepInternal() {
        assert(Thread::getId()==0);
        assert(instance_->bUseDefaultSchedule_ == true);
        // Allow multiple calls, but ignore if schedule is active.
        if (isActive_) {
            assert(instance_ == this);
            return;
        }
        // Cannot swap out an active schedule (call wait first)
        assert(instance_->isActive_ == false);
        instance_ = this;
        isActive_ = true;
        for (auto &task: tasks_) {
            task.functor_ = nullptr;
            task.functorBeginBarrier_.close();
            task.functorEndBarrier_.close(task.getNumThreads());
        }

        // Increment counter only
        stepCounter_.fetch_add(1, std::memory_order_release);
    }
    //! Wait on all tasks to finish
    void waitInternal() {
        assert(Thread::getId()==0);
        assert(this == instance_);
        // Allow multiple calls, but ignore if schedule is not active.
        if (!isActive_) {
            return;
        }
        threads_[0].processQueue(); //Before waiting the OS thread executes its queue
        for(unsigned int i=1;i<tasks_.size();i++) {
            tasks_[i].functorEndBarrier_.wait();
        }
        isActive_ = false;
        instance_ = defaultInstance_;
    }
    std::deque<Task>  tasks_;  //It is essential this isn't a vector (doesn't get moved when resizing). Is this ok to be a list (linear) or does it need to be a tree? A serial task isn't before a loop. It is both before and after.
    std::map<std::string,int> taskLabels_;
    std::vector< std::vector<SubTask *> > threadSubTasks_;
    bool bUseDefaultSchedule_ = false;
    // "Active" means schedule is between nextStep and wait calls.
    bool isActive_ = false;
    static std::deque<Thread> threads_;
    static std::atomic<int> stepCounter_;
    static STS *defaultInstance_;
    static std::map<std::string, STS *> stsInstances_;
    static STS *instance_;
};

#endif // STS_STS_H
