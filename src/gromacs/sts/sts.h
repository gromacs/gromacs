#ifndef STS_STS_H
#define STS_STS_H

#include <cassert>

#include <atomic>
#include <deque>
#include <map>
#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <memory>

#include "range.h"
#include "reduce.h"
#include "task.h"
#include "thread.h"

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
        // Users cannot create an instance with the default schedule but can
        // use this default instance as a quick way to parallelize loops.
        defaultInstance_ = new STS("default");
        for (int id = 0; id < numThreads; id++) {
            defaultInstance_->assign("default", id, {{id, numThreads},
                                                    {id + 1, numThreads}});
        }
        defaultInstance_->bUseDefaultSchedule_ = true;
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
    //! Notify threads to start computing the next step
    void nextStep() {
        assert(Thread::getId()==0);
        assert(instance_->bUseDefaultSchedule_ == true);
        instance_ = this;
        for (auto &task: tasks_) {
            task.functor_ = nullptr;
            task.functorBeginBarrier_.close();
            task.functorEndBarrier_.close(task.getNumThreads());
        }
        stepCounter_.fetch_add(1, std::memory_order_release);
    }
    /*! \brief
     * Run an asynchronous function
     *
     * \param[in] label     The task label (needs to match assign())
     * \param[in] function  The function (or lambda) to execute
     */
    template<typename F>
    void run(std::string label, F function) {
        assert(this == instance_);
        if (!isTaskAssigned(label) || bUseDefaultSchedule_) {
            function();
        } else {
            tasks_[getTaskId(label)].functor_ = new BasicTaskFunctor<F>(function);
            tasks_[getTaskId(label)].functorBeginBarrier_.open();
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
        assert(this == instance_);
        int taskId = -1;
        if (bUseDefaultSchedule_) {
            taskId = 0;
            nextStep(); //Default schedule has only a single step and the user doesn't need to call nextStep
            assert(getTaskId("default")==taskId);
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
        int nextSubTaskId = thread.getCurrentSubTaskId() + 1;
        assert(getSubTask(Thread::getId(), nextSubTaskId)->getTaskId() == taskId);
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
            wait();
            // Calling wait sets the current instance to default, which we
            // don't want to do internally but only at user level.
            instance_ = this;
        }
    }
    //! Automatically compute new schedule based on previous step timing
    void reschedule() {
        // not yet available
    }
    //! Wait on all tasks to finish
    void wait() {
        assert(this == instance_);
        threads_[0].processQueue(); //Before waiting the OS thread executes its queue
        for(unsigned int i=1;i<tasks_.size();i++) {
            tasks_[i].functorEndBarrier_.wait();
        }
        if (bSTSDebug_) {
            for (const auto &t : tasks_) {
                for (const auto &st : t.subtasks_) {
                    auto wtime = std::chrono::duration_cast<std::chrono::microseconds>(st->waitTime_).count();
                    auto rtime = std::chrono::duration_cast<std::chrono::microseconds>(st->runTime_).count();
                }
                if (t.subtasks_.size() > 1) {
                    auto ltwtime = std::chrono::duration_cast<std::chrono::microseconds>(t.waitTime_).count();
                }
            }
        }
        instance_ = defaultInstance_;
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
        return TaskReduction<T>(taskName, init, numThreads);
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
    bool isTaskAssigned(std::string label) {
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
    /*! \brief
     * Returns task label for task ID
     *
     * \param[in] id   task Id
     * \returns        task label
     */
    std::string getTaskLabel(int id) const {
        for (auto it: taskLabels_) {
            if (it.second == id) return it.first;
        }
        throw std::invalid_argument("Invalid task Id: "+id);
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
    std::deque<Task>  tasks_;  //It is essential this isn't a vector (doesn't get moved when resizing). Is this ok to be a list (linear) or does it need to be a tree? A serial task isn't before a loop. It is both before and after.
    std::map<std::string,int> taskLabels_;
    std::vector< std::vector<SubTask *> > threadSubTasks_;
    bool bUseDefaultSchedule_ = false;
    bool bSTSDebug_ = true;
    static std::deque<Thread> threads_;
    static std::atomic<int> stepCounter_;
    static STS *defaultInstance_;
    static std::map<std::string, STS *> stsInstances_;
    static STS *instance_;
};

#endif // STS_STS_H
