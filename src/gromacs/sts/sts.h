#ifndef STS_STS_H
#define STS_STS_H

#include <cassert>

#include <algorithm>
#include <deque>
#include <iostream>
#include <map>
#include <memory>
#include <stack>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <utility>
#include <vector>

#include <atomic>
#include <chrono>

#include "range.h"
#include "reduce.h"
#include "task.h"
#include "thread.h"

/*! \internal \brief
 * Static Thread Scheduler
 */
class STS {
public:
    /*! \brief Schedule Id
     *
     * Can be used to retrieve schedules with getInstance() rather than storing
     * them in the application code.
     */
    const std::string id;
    /*! \brief
     * Startup STS and set the number of threads.
     * No STS functions should be called before startup.
     */
    static void startup(size_t numThreads) {
        assert(numThreads > 0);
        if (threads_.size() > 0) {
            // Ignore multiple calls but check that thread count is the same
            assert(numThreads == threads_.size());
            return;
        }

        // Barrier must be initialized before creating threads
        // First time the barrier is used, each non-main thread enters it twice
        stepCounterBarrier_.close(2*(numThreads-1));

        for (size_t id = 0; id < numThreads; id++) {
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
        for (int i=1;i<getNumThreads();i++) {
            threads_[i].join();
        }
        delete defaultInstance_;
    }
    /*! \brief
     * Constructs a new STS schedule
     *
     * \param[in] name  optional name for schedule
     */
    STS(std::string name = "") :id(name), bUseDefaultSchedule_(false), isActive_(false) {
        int n = getNumThreads();
        assert(n > 0);
        threadSubTasks_.resize(n);
        threadCallStacks_.resize(n);
        if (!id.empty()) {
            stsInstances_[id] = this;
        }
    }
    ~STS() {
        // Tasks own their subtasks, so do not delete SubTasks in threadSubTasks_
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
     * This is the most generic assign function, which does the actual work.
     * Most of the time, applications should use one of the more specific
     * assign functions.
     *
     * Note that this function can only assign one thread at a time, whereas
     * some other assign functions can assign a range of threads.
     *
     * If a range for a loop task is specified, only that section of the loop is assigned.
     * In that case it is important to assign the remaining loop out of [0,1] also to
     * some other thread. It is valid to assign multiple parts of a loop to the same thread.
     * The order of assign calls specifies in which order the thread executes the tasks.
     *
     * \param[in] label    Label of the task. Needs to match the run()/parallel_for() label
     * \param[in] ttype    Task type
     * \param[in] threadId Id of thread assigned the work
     * \param[in] range    Assigned range for a loop task. Ignored for a basic task.
     */
    enum class TaskType {BASIC, LOOP, MULTILOOP};
    void assign(std::string label, TaskType ttype, int threadId, Range<Ratio> range) {
        int id = setTask(label, ttype);
        assert(range.start>=0 && range.end<=1);
        SubTask* t = new SubTask(tasks_[id].get(), range);
        threadSubTasks_[threadId].push_back(t);
        tasks_[id]->pushSubtask(threadId, t);
    }
    /*! \brief
     * Assign a basic task to a single thread
     *
     * \param[in] label    Label of the task
     * \param[in] threadId Id of thread assigned the work
     */
    void assign_run(std::string label, int threadId) {
        assign(label, TaskType::BASIC, threadId, Range<Ratio>(1));
    }
    /*! \brief
     * Assign a basic task to a single thread along with a set of helper threads
     * to execute contained loops.
     *
     * Note that helper threads wait until the basic task completes. Loops should
     * be assigned explicitly if more fine-grained control is needed.
     *
     * \param[in] label         Label of the task
     * \param[in] threadId      Id of thread assigned the work
     * \param[in] helperThreads vector of threads that execute contained loops
     */
    void assign_run(std::string label, int threadId, const std::vector<int> &helperThreads) {
        assign_run(label, threadId);
        int nthreads = helperThreads.size();
        int numer = 0;
        int denom = nthreads;
        std::string loopTaskLabel = label + "_multiloop";
        // Main thread of a basic task is always a helper thread too. So add it
        // as a helper if not already listed.
        if (std::find(helperThreads.begin(), helperThreads.end(), threadId) == helperThreads.end()) {
            denom++;
            numer++;
            assign(loopTaskLabel, TaskType::MULTILOOP, threadId, {1, denom});
        }
        for (int i=0; i<nthreads; i++, numer++) {
            assign(loopTaskLabel, TaskType::MULTILOOP, helperThreads[i], {{numer, denom},{numer+1, denom}});
        }

        // Link the basic task with its loop task
        int parentId = getTaskId(label);
        int childId  = getTaskId(loopTaskLabel);
        Task* pTask = tasks_[parentId].get();
        Task* cTask = tasks_[childId].get();
        assert(std::type_index(typeid(*pTask))==std::type_index(typeid(BasicTask)));
        assert(std::type_index(typeid(*cTask))==std::type_index(typeid(MultiLoopTask)));
        BasicTask*     parentTask = static_cast<BasicTask*>(pTask);
        MultiLoopTask* childTask  = static_cast<MultiLoopTask*>(cTask);
        parentTask->setMultiLoop(childTask);
    }
    /*! \brief
     * Assign a loop task to a single thread
     *
     * \param[in] label    Label of the task
     * \param[in] threadId Id of thread assigned the work
     * \param[in] range    range of loop executed by this thread
     */
    void assign_loop(std::string label, int threadId, const Range<Ratio> range) {
        assign(label, TaskType::LOOP, threadId, range);
    }
    /*! \brief
     * Assign a loop task to a vector of threads. This will split the loop
     * evenly among the given threads.
     *
     * \param[in] label     Label of the task
     * \param[in] threadIds Id of threads assigned the work
     */
    void assign_loop(std::string label, const std::vector<int> &threadIds) {
        int nthreads = threadIds.size();
        for (int i=0; i<nthreads; i++) {
            assign_loop(label, threadIds[i], {{i, nthreads},{i+1, nthreads}});
        }
    }
    void setNormalPriority(std::string label) {
        tasks_[getTaskId(label)]->setPriority(Task::NORMAL);
    }
    template<typename ... Types>
    void setNormalPriority(std::string label, Types ... rest) {
        tasks_[getTaskId(label)]->setPriority(Task::NORMAL);
        setNormalPriority(rest...);
    }
    void setHighPriority(std::string label) {
        tasks_[getTaskId(label)]->setPriority(Task::HIGH);
    }
    template<typename ... Types>
    void setHighPriority(std::string label, Types ... rest) {
        tasks_[getTaskId(label)]->setPriority(Task::HIGH);
        setHighPriority(rest...);
    }
    //! \brief Clear all assignments
    void clearAssignments() {
        for (auto &taskList : threadSubTasks_) {
            taskList.clear();
        }
        for (std::unique_ptr<Task> &task : tasks_) {
            task->clearSubtasks();
        }
    }
    /*! \brief
     * Set this schedule to use the default schedule
     *
     * All previous assignments are cleared. Loops are divided evenly among all threads and
     * non-loop tasks simply run on the invoking thread,
     */
    void setDefaultSchedule() {
        bUseDefaultSchedule_ = true;
        clearAssignments();
        int numThreads = getNumThreads();
        for (int id = 0; id < numThreads; id++) {
            assign_loop("default", id, {{id, numThreads}, {id + 1, numThreads}});
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
            tasks_[getTaskId(label)]->setFunctor(new BasicTaskFunctor<F>(function));
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
            // Task object needed to do reductions
            assert(red == nullptr);
            for (int i=start; i<end; i++) {
                body(i);
            }
            return;
        } else {
            taskId = getTaskId(label);
        }
        Task* task = tasks_[taskId].get();
        assert(task != nullptr);
        // Must be a loop task
        assert(std::type_index(typeid(*task))!=std::type_index(typeid(BasicTask)));
        task->setReduction(red);
        task->setFunctor(new LoopTaskFunctor<F>(body, {start, end}));
        runNestedLoop(task);
        if (red != nullptr) {
            red->reduce();
        }
        // User does not need to call wait for default scheduling
        if (bUseDefaultSchedule_) {
            waitInternal();
        }
    }
    /*! \brief
     * Skip the given task.
     *
     * This is useful when an assigned task should not run under certain conditions.
     *
     * \param[in] label  task label
     */
    void skipTask(std::string label) {
        if (isTaskAssigned(label)) {
            int taskId = getTaskId(label);
            Task* task = tasks_[taskId].get();
            assert(task != nullptr);
            task->markDone();
        }
    }
    //! Automatically compute new schedule based on previous step timing
    void reschedule() {
        // not yet available
    }
    //! Wait for specific task to finish
    void waitForTask(std::string label) {
        assert(isTaskAssigned(label));
        int t = getTaskId(label);
        tasks_[t]->wait();
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
     * \param[in] id  STS instance Id
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
     * WARNING: meant only for internal use. Applications should use
     * "getInstance" for better error checking and clarity when using
     * multiple STS instances.
     *
     * \returns current STS instance
     */
    static STS *getCurrentInstance() {
        return instance_;
    }
    /*! \brief
     * Get number of subtasks assigned to a thread
     *
     * \param[in] threadId  Thread Id
     */
    int getNumSubTasks(int threadId) const {
        return threadSubTasks_[threadId].size();
    }
    /*! \brief
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
    /*! \brief
     * Get number of threads for a given task
     *
     * \param[in] label  task label
     * \return number of threads
     */
    int getTaskNumThreads(std::string label) {
        assert(isTaskAssigned(label));
        int taskId = getTaskId(label);
        return tasks_[taskId]->getNumThreads();
    }
    /*! \brief
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
    /*! \brief
     * Load atomic step counter
     *
     * \returns step counter
     */
    static int loadStepCounter() { return stepCounter_.load(std::memory_order_acquire); }
    /*! \brief
     * Wait on atomic step counter to change
     *
     * param[in] c   last step processed by thread
     */
    static int waitOnStepCounter(int c) {
        stepCounterBarrier_.markArrival();
        return wait_until_not(stepCounter_, c);
    }

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
    bool runNextSubTask() {
        int tid = Thread::getId();
        // Should only be called to start a new call stack
        assert(threadCallStacks_[tid].empty());
        for (int stid=0; stid < getNumSubTasks(tid); stid++) {
            if (!threadSubTasks_[tid][stid]->isDone()) {
                runSubTask(stid);
                return true;
            }
        }
        return false;
    }
    /*! \brief
     * Yield a running task and run the next high priority task
     */
    void yield() {
        int tid = Thread::getId();
        const std::stack<int>& cstack = threadCallStacks_[Thread::getId()];
        int stid = 0;
        if (!cstack.empty()) {
            stid = cstack.top()+1;
        }
        for (; stid < getNumSubTasks(tid); stid++) {
            SubTask* subtask = threadSubTasks_[tid][stid];
            // TODO: Decide on the exact policy for what tasks will be considered
            if (!subtask->isDone() && subtask->getTask()->getPriority() ==
                    Task::HIGH && subtask->isReady()) {
                runSubTask(stid);
            }
        }
    }
private:
    const Task* getCurrentTask() {
        int tid = Thread::getId();
        if (threadCallStacks_[tid].empty()) {
            return nullptr;
        }
        else {
            int stid = threadCallStacks_[tid].top();
            return threadSubTasks_[tid][stid]->getTask();
        }
    }
    bool isTaskAssigned(std::string label) const {
        return (taskLabels_.find(label) != taskLabels_.end());
    }
    /* \brief
     * Get task id for task label
     *
     * Undefined behavior if task doesn't exist.
     *
     * \param[in] task label
     * \return task id
     */
    int getTaskId(std::string label) {
        auto it = taskLabels_.find(label);
        assert(it != taskLabels_.end());
        return it->second;
    }
    /* \brief
     * Sets values for a task, creating a new task object if it doesn't exist.
     *
     * Creating tasks isn't thread safe, so this should only be called when
     * doing thread assignments between schedule runs.
     * \param[in] label for task
     * \return task id
     */
    int setTask(std::string label, TaskType ttype) {
        // TODO: Add asserts for schedule state (using isActive_ variable perhaps)
        assert(Thread::getId() == 0);
        auto it = taskLabels_.find(label);
        if (it != taskLabels_.end()) {
            return it->second;
        }
        else {
            unsigned int v = taskLabels_.size();
            assert(v==tasks_.size());
            switch (ttype) {
            case TaskType::BASIC:
                tasks_.emplace_back(std::unique_ptr<Task>(new BasicTask()));
                break;
            case TaskType::LOOP:
                tasks_.emplace_back(std::unique_ptr<Task>(new LoopTask()));
                break;
            case TaskType::MULTILOOP:
                tasks_.emplace_back(std::unique_ptr<Task>(new MultiLoopTask()));
                break;
            }
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
        const Task* t = getCurrentTask();
        // Must be a loop task
        assert(std::type_index(typeid(*t))!=std::type_index(typeid(const BasicTask)));
        // TODO: How to handle these user errors - calling collect outside of a
        // task or for a task without a reduction?
        if (t == nullptr) {
            return;
        }
        TaskReduction<T>* tr = (static_cast<TaskReduction<T> *>(t->getReduction()));
        if (tr == nullptr) {
            return;
        }
        tr->collect(a, ttid);
    }
    /*! \brief
     * Run the given subtask
     *
     * This function handles the low-level operations of running a subtask.
     * This includes bookkeeping, setting flags, and launching the run.
     * Callers are responsible for selecting a suitable subtask to run.
     *
     * \param[in] stid  Id of the subtask to run
     */
    void runSubTask(int stid) {
        int tid = Thread::getId();
        std::stack<int>& cstack = threadCallStacks_[tid];
        cstack.push(stid);
        bool isDone = threadSubTasks_[tid][stid]->run();
        cstack.pop();
        threadSubTasks_[tid][stid]->setDone(isDone);
    }
    /*! \brief
     * Run the given loop task inside the currently running task.
     *
     * The caller is the main thread for the loop, not a helper, and thus the
     * task must be the next unfinished task on the caller's queue. Also, the
     * caller waits on the task to finish at the end of the loop.
     *
     * \param[in] task  The loop task to be executed.
     */
    void runNestedLoop(Task* task) {
        int tid = Thread::getId();
        const std::stack<int>& cstack = threadCallStacks_[Thread::getId()];
        int stid = 0;
        if (!cstack.empty()) {
            stid = cstack.top()+1;
        }
        for (; stid < getNumSubTasks(tid); stid++) {
            SubTask* subtask = threadSubTasks_[tid][stid];
            if (!subtask->isDone()) {
                // Loop task must be the thread's next unfinished task
                assert(subtask->getTask() == task);
                runSubTask(stid);
                task->wait();
                return;
            }
        }
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
        for (std::unique_ptr<Task> &task: tasks_) {
            task->restart();
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
            tasks_[i]->wait();
        }

        // Wait for all threads to complete step before changing any internal state
        stepCounterBarrier_.wait();
        stepCounterBarrier_.close(getNumThreads()-1);

        isActive_ = false;
        instance_ = defaultInstance_;
    }
    // TODO: Consider using a tree of tasks if scheduling becomes unwieldy.
    // Nesting of loops inside basic tasks already makes the logic somewhat hard
    // to follow with this simple list.
    std::vector<std::unique_ptr<Task>>  tasks_;
    std::map<std::string,int> taskLabels_;
    std::vector< std::vector<SubTask *> > threadSubTasks_;
    // Call stack of running subtasks for each thread
    std::vector<std::stack<int>> threadCallStacks_;
    bool bUseDefaultSchedule_;
    // "Active" means schedule is between nextStep and wait calls.
    bool isActive_;
    // Cannot be a vector because Task moving is not allowed (occurs on resizing)
    static std::deque<Thread> threads_;
    static std::atomic<int> stepCounter_;
    static OMBarrier stepCounterBarrier_;
    static STS *defaultInstance_;
    static std::map<std::string, STS *> stsInstances_;
    static STS *instance_;
};

#endif // STS_STS_H
