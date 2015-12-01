#ifndef STS_THREAD_H
#define STS_THREAD_H

#include <cassert>

#include <atomic>
#include <chrono>
#include <deque>
#include <memory>
#include <thread>

#include "range.h"

#if (__GNUC__ == 4 && __GNUC_MINOR__ <= 7) || (defined __ICC && __ICC <= 1400)
#define thread_local __thread
#endif
#if _MSC_VER == 1800
#define thread_local __declspec(thread)
#endif

using sts_clock = std::chrono::steady_clock;

/*! \brief
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
    SubTask(int taskId, Range<Ratio> range) : taskId_(taskId), range_(range) {
        setDone(false);
    }
    //! Wait on sub-task to complete
    void wait() const
    {
        while(!(done_.load(std::memory_order_acquire))); //TODO: add gmx_pause
    }
    /*! \brief
     * Set done status.
     *
     * \param[in] done   value
     */
    void setDone(bool done)
    {
        done_.store(done, std::memory_order_release);
    }
    int getTaskId() const { return taskId_; }               //!< get Task Id
    const Range<Ratio>& getRange() const { return range_; } //!< get Range

    sts_clock::duration waitTime_; /**< Time spent until task was ready */
    sts_clock::duration runTime_;  /**< Time spent executing sub-task  */
private:
    int taskId_;                   /**< The ID of the task this is a part of */
    Range<Ratio> range_;           /**< Range (out of [0,1]) of loop part */
    std::atomic_bool done_;        /**< Sub-task is done */
};

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
    Thread(int id) {
        if (id!=0) {
            thread_.reset(new std::thread([=](){id_=id; doWork();}));
        }
    }
    Thread(Thread&&) = delete; //nextSubtaskId_ access would not be thread-safe
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
    /*! \brief
     * Add a subtask to the end of the queue
     *
     * \param[in] taskID  The ID of the task to add
     * \param[in] range   The range of the loop to be done by this thread
     *                    Ignored for basic tasks.
     * \returns           Pointer to sub-task added
    */
    SubTask const* addSubtask(int taskId, Range<Ratio> range) {
        taskQueue_.emplace_back(taskId, range);
        return &(taskQueue_.back());
    }
    /*! \brief
     * Clear all sub-tasks in the queue.
     *
     * Is called by clearAssignments to prepare for rescheduling tasks. */
    void clearSubtasks() { taskQueue_.clear(); }
    /*! \brief
     * Return sub-task next in queue
     *
     * \returns      Pointer to sub-task
     */
    SubTask const* getNextSubtask() { return &(taskQueue_[nextSubtaskId_]); }
    /*! \brief
     * Resets queue in prepartion for the next step
     *
     * Does not delete queue. @see clearSubtasks
     */
    void resetTaskQueue() {
        if (!thread_) nextSubtaskId_ = 0;
        for (auto &s: taskQueue_) s.setDone(false);
    }
    //! Wait for thread to finish
    void join() { if (thread_) thread_->join(); }
    /*! \brief
     * Return thread Id
     *
     * Note: This is a static method. Id returned depends on the thread executing
     * not the object. Only use Thread::getID() not t.getId() to avoid confusion.
     *
     * \returns thread id
     */
    static int getId() { return id_; }
    //! Wait for thread to finish all tasks in queue (=step)
    void wait() {
        for(auto &task: taskQueue_) task.wait(); //If we have a task-graph it is sufficient to wait on last parent task. Without graph we need to make sure all are done.
    }
private:
    void doWork(); //function executed by worker threads

    std::deque<SubTask> taskQueue_;
    unsigned int nextSubtaskId_ = 0;
    std::unique_ptr<std::thread> thread_;
    static thread_local int id_;
};

#endif // STS_THREAD_H
