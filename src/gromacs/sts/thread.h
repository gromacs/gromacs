#ifndef STS_THREAD_H
#define STS_THREAD_H

#include <cassert>

#include <atomic>
#include <deque>
#include <memory>
#include <thread>

#include "range.h"
#include "task.h"

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
    Thread(int id) {
        if (id!=0) {
            thread_.reset(new std::thread([=](){id_=id; doWork();}));
        }
        //setAffinity(id);
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
     * \returns thread Id
     */
    static int getId() { return id_; }
    int getCurrentTaskId() {
        if (nextSubtaskId_ < 1) {
            return -1;
        }
        return taskQueue_[nextSubtaskId_ - 1].getTaskId();
    }
    //! Wait for thread to finish all tasks in queue (=step)
    void wait() {
        for(auto &task: taskQueue_) task.wait(); //If we have a task-graph it is sufficient to wait on last parent task. Without graph we need to make sure all are done.
    }
private:
    void doWork(); //function executed by worker threads
    void setAffinity(int coreId);

    std::deque<SubTask> taskQueue_;
    unsigned int nextSubtaskId_ = 0;
    std::unique_ptr<std::thread> thread_;
    static thread_local int id_;
};

#endif // STS_THREAD_H
