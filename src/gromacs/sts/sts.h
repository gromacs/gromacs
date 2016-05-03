#ifndef STS_STS_H
#define STS_STS_H

#include <cassert>

#include <atomic>
#include <deque>
#include <map>
#include <string>
#include <chrono>
#include <iostream>
#include <memory>

#include "range.h"
#include "reduce.h"
#include "task.h"
#include "thread.h"

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
    /*! \brief
     * Constructor
     */
    STS() {
        stepCounter_.store(0, std::memory_order_release);
        threads_.emplace_back(0);
    }
    ~STS() {
        //-1 notifies threads to finish
        stepCounter_.store(-1, std::memory_order_release);
        for(unsigned int i=1;i<threads_.size();i++) {
            threads_[i].join();
        }
    }
    /*! \brief
     * Set number of threads in the pool
     *
     * \param[in] n number of threads to use (including OS thread)
     */
    void setNumThreads(int n) {
        for (int id = threads_.size(); id < n; id++) {
            threads_.emplace_back(id); //create threads
        }
        for (int id = threads_.size(); id > n; id--) {
            threads_.pop_back();
        }
        if (bUseDefaultSchedule_) {
            clearAssignments();
            for (int id = 0; id < n; id++) {
                assign("default", id, {{id,     n},
                        {id + 1, n}});
            }
            bUseDefaultSchedule_ = true;
        }
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
    void assign(std::string label, int threadId, Range<Ratio> range) {
        int id = getTaskId(label);
        assert(range.start>=0 && range.end<=1);
        SubTask const* subtask = threads_.at(threadId).addSubtask(id, range);
        tasks_[id].pushSubtask(threadId, subtask);
        bUseDefaultSchedule_ = false;
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
    //! Clear all assignments
    void clearAssignments() {
        for (auto &thread : threads_) {
            thread.clearSubtasks();
        }
        for (auto &task : tasks_) {
            task.clearSubtasks();
        }
    }
    //! Notify threads to start computing the next step
    void nextStep() {
        assert(Thread::getId()==0);
        for (auto &task: tasks_) {
            task.functor_.reset(nullptr);
        }
        for (int i=0; i<threads_.size(); i++) {
            threads_[i].resetTaskQueue();
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
        if (bUseDefaultSchedule_) {
            function();
        } else {
            tasks_[getTaskId(label)].functor_.reset(new BasicTaskFunctor<F>(function));
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
        int taskId = 0;
        if (bUseDefaultSchedule_) {
            nextStep(); //Default schedule has only a single step and the user doesn't need to call nextStep
            assert(getTaskId("default")==taskId);
        } else {
            taskId = getTaskId(label);
        }
        auto &task = tasks_[taskId];
        task.reduction_ = red;
        task.functor_.reset(new LoopTaskFunctor<F>(body, {start, end}));
        auto &thread = threads_[Thread::getId()];
        //Calling processTask implies that the thread calling parallel_for participates in the loop and executes it next in queue
        assert(thread.getNextSubtask()->getTaskId()==taskId);
        thread.processTask();
        auto startWaitTime = sts_clock::now();
        for(auto s: task.subtasks_) {
            s->wait();
        }
        task.waitTime_ = sts_clock::now() - startWaitTime;

        // TODO: A smarter reduction would take place before the above wait.
        if (task.reduction_ != nullptr) {
            auto startReductionTime = sts_clock::now();
            static_cast< TaskReduction<T> *>(task.reduction_)->reduce();
            task.reductionTime_ = sts_clock::now() - startReductionTime;
        }
    }
    //! Automatically compute new schedule based on previous step timing
    void reschedule() {
        // not yet available
    }
    //! Wait on all tasks to finish
    void wait() {
        if (!bUseDefaultSchedule_) {
            threads_[0].processQueue(); //Before waiting the OS thread executes its queue
            for(unsigned int i=1;i<threads_.size();i++) {
                threads_[i].wait();
            }
            if (bSTSDebug_) {
                std::cerr << "Times for step " << loadStepCounter() << std::endl;
                for (const auto &t : tasks_) {
                    for (const auto &st : t.subtasks_) {
                        auto wtime = std::chrono::duration_cast<std::chrono::microseconds>(st->waitTime_).count();
                        auto rtime = std::chrono::duration_cast<std::chrono::microseconds>(st->runTime_).count();
                        std::cerr << getTaskLabel(st->getTaskId()) << " " << wtime << " " << rtime << std::endl;
                    }
                    if (t.subtasks_.size() > 1) {
                        auto ltwtime = std::chrono::duration_cast<std::chrono::microseconds>(t.waitTime_).count();
                        std::cerr << "Wait for task to complete " << ltwtime << std::endl;
                    }
                }
            }
        }
    }
    /*! \brief
     * Returns the STS instance
     *
     * \returns STS instance
     */
    static STS *getInstance() { return instance_.get(); }
    /*! \brief
     * Returns the task functor for a given task Id
     *
     * Waits on functor to be ready if the corresponding run()/parallel_for() hasn't been executed yet.
     *
     * \param[in] task Id
     * \returns task functor
     */
    ITaskFunctor *getTaskFunctor(int taskId) {
        auto &func = tasks_[taskId].functor_;
        return func.wait();
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
        // TODO: Handle case where label is not a valid task.
        // Currently, it will insert a new task!
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
    int loadStepCounter() { return stepCounter_.load(std::memory_order_acquire); }
    /* \brief
     * Wait on atomic step counter to change
     *
     * param[in] c   last step processed by thread
     */
    int waitOnStepCounter(int c) {return wait_until_not(stepCounter_, c);}
private:
    // Helper function for operations that need the current task
    const Task *getCurrentTask() {
        int threadId = Thread::getId();
        int taskId = threads_[threadId].getCurrentTaskId();
        if (taskId == -1) {
            return nullptr;
        }
        return &tasks_[taskId];
    }
    //Creates new ID for unknown label.
    //Creating IDs isn't thread safe. OK because assignments and run/parallel_for (if run without pre-assignment) are executed by master thread while other threads wait on nextStep.
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

    std::deque<Task>  tasks_;  //It is essential this isn't a vector (doesn't get moved when resizing). Is this ok to be a list (linear) or does it need to be a tree? A serial task isn't before a loop. It is both before and after.
    std::map<std::string,int> taskLabels_;
    std::deque<Thread> threads_;
    static std::unique_ptr<STS> instance_;
    bool bUseDefaultSchedule_ = true;
    bool bSTSDebug_ = true;
    std::atomic<int> stepCounter_;
};

#endif // STS_STS_H
