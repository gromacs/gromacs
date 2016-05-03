#ifndef STS_TASK_H
#define STS_TASK_H

#include <atomic>
#include <chrono>
#include <map>
#include <vector>

#include "range.h"
#include "wait.h"

using sts_clock = std::chrono::steady_clock;

//! Interface of the executable function of a task
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

//! Loop task functor
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

//! Basic (non-loop) task functor
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

//! Atomic scoped_ptr like smart pointer
template<class T>
class AtomicPtr {
public:
    AtomicPtr() { ptr_.store(nullptr, std::memory_order_release); }
    ~AtomicPtr() {
        T* p = get();
        if(p) delete p;
    }
    /*! \brief
     * Deletes the previous object and stores new pointer
     *
     * \param[in] t   pointer to store
     */
    void reset(T* t) {
        T* o = ptr_.exchange(t, std::memory_order_release);
        if(o) delete o;
    }
    /*! \brief
     * Returns the stored pointer
     *
     * \returns pointer
     */
    T* get() { return ptr_.load(std::memory_order_consume); }
    /*! \brief
     * Wait for pointer to be set (non-null)
     *
     * \returns pointer
     */
    T* wait() { return wait_until_not(ptr_, static_cast<T*>(nullptr)); }
private:
    std::atomic<T*> ptr_;
};

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
        wait_until(done_, true);
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
 * A task to be executed
 *
 * Can either be a function or loop. Depending on the schedule is
 * executed synchronous or asynchronous. Functions are always
 * executed by a single thread. Loops are executed, depending on
 * the schedule, in serial or in parallel.
 */
struct Task {
    AtomicPtr<ITaskFunctor> functor_;      //!< The function/loop to execute
    void *reduction_;
    //! All subtasks of this task. One for each section of a loop. One for a basic task.
    std::vector<SubTask const*> subtasks_;
    //!< The waiting time in the implied barrier at the end of a loop. Zero for basic task.
    sts_clock::duration waitTime_;
    sts_clock::duration reductionTime_;

    Task() :reduction_(nullptr), waitTime_(0), reductionTime_(0), numThreads_(0) {
        functor_.reset(nullptr);
    }
    void pushSubtask(int threadId, SubTask const* t) {
        subtasks_.push_back(t);
        if (threadTaskIds_.find(threadId) == threadTaskIds_.end()) {
            threadTaskIds_[threadId] = numThreads_;
            numThreads_++;
        }
    }
    void clearSubtasks() {
        subtasks_.clear();
        threadTaskIds_.clear();
        numThreads_ = 0;
    }
    int getNumThreads() const {
        return numThreads_;
    }
    int getThreadId(int threadId) const {
        auto id = threadTaskIds_.find(threadId);
        if (id == threadTaskIds_.end()) {
            return -1;
        }
        return (*id).first;
    }
private:
    int numThreads_;
    //! Map global thread id to an id only for this task (task ids are consecutive starting from 0)
    std::map<int, int> threadTaskIds_;
};

#endif // STS_TASK_H
