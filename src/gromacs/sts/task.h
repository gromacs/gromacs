#ifndef STS_TASK_H
#define STS_TASK_H

#include <map>
#include <set>
#include <vector>

#include <chrono>

#include "barrier.h"
#include "range.h"
#include "thread.h"

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

class Task;

/*! \internal \brief
 * The portion of a task done by one thread
 *
 * Contains run() method to execute subtask directly, along with needed
 * information (task and range).
 */
class SubTask {
public:
    /*! \brief
     * Constructor
     *
     * \param[in] task     The task this is part of.
     * \param[in] range    Out of a possible range from 0 to 1, the section in
                           this part. Ignored for basic tasks.
     */
    SubTask(Task *task, Range<Ratio> range) :range_(range), task_(task),
    isDone_(false) {}
    /*! \brief
     * Run the subtask
     */
    bool run();
    const Task *getTask() const;
    bool isDone() const;
    void setDone(bool isDone);
    bool isReady() const;

    const Range<Ratio> range_;     /**< Range (out of [0,1]) of loop part */
private:
    Task *task_;             /**< Reference to main task */
    bool isDone_;
};

/*! \internal \brief
 * Base task class that handles the storage and bookkeeping of subtasks,
 * common operations needed by all subclasses. It also stores reduction
 * functions and task priority.
 *
 * Each subtask is assigned to a single thread, and this class assigns a
 * Task-specific thread id to all participating threads. Subclasses store
 * the actual function and actually execute the task.
 *
 * Note that for non-loop tasks, only a single subtask and thread are
 * needed, and thus much of this infrastructure is unused.
 */
class Task {
public:
    enum Priority {NORMAL, HIGH};
    Task() :reduction_(nullptr), numThreads_(0), priority_(NORMAL) {}
    /*! \brief
     * Add a new subtask for this task
     *
     * Note: Task takes ownership of passed subtask.
     *
     * \param[in] threadId  thread to which this subtask is assigned
     * \param[in] t         subtask
     */
    void pushSubtask(int threadId, SubTask* t) {
        subtasks_.emplace_back(t);
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
        return (*id).second;
    }
    /*! \brief
     * Restart the task. Must be called on each task prior to starting work,
     * normally called for all tasks in an STS schedule at the beginning of
     * each step.
     */
    void restart() {
        for (std::unique_ptr<SubTask> &st : subtasks_) {
            st->setDone(false);
        }
        init();
    }
    /*! \brief
     * Restart the task. Must be called on each task prior to starting work,
     * normally called for all tasks in an STS schedule at the beginning of
     * each step.
     */
    void markDone() {
        for (std::unique_ptr<SubTask> &st : subtasks_) {
            st->setDone(true);
        }
    }
    //! \brief Get task priority
    Priority getPriority() const {
        return priority_;
    }
    //! \brief Set task priority
    void setPriority(Priority p) {
        priority_ = p;
    }
    //! \brief Get reduction function
    void* getReduction() const {
        return reduction_;
    }
    /*! \brief
     * Store a new reduction function
     *
     * This class only stores the reduction. Clients and subclasses are
     * responsible for making use of the reduction.
     *
     * \param[in] r reduction function
     */
    void setReduction(void* r) {
        reduction_ = r;
    }
    //! \brief Set the functor (work) to be done
    virtual void setFunctor(ITaskFunctor*) = 0;
    //! \brief Return whether task is ready to run
    virtual bool isReady() const = 0;
    //! \brief Run the functor in the specified range
    virtual bool run(Range<Ratio>) = 0;
    //! \brief Wait for the functor to complete
    virtual void wait() = 0;
    virtual ~Task() {}
protected:
    void*    reduction_;
private:
    /*! \brief
     * Initialize the task. This function is called by "restart" and should do
     * all necessary steps to ready the task for doing work (initializing
     * barriers, counters, variables, etc.).
     */
    virtual void  init() = 0;
    //! All subtasks of this task. One for each section of a loop. One for a basic task.
    std::vector<std::unique_ptr<SubTask>> subtasks_; //!< Subtasks to be executed by a single thread
    int numThreads_;
    //! Map STS thread id to an id only for this task (task ids are consecutive starting from 0)
    std::map<int, int> threadTaskIds_;
    Priority priority_;
};

/*! \brief
 * Class for loop tasks
 */
class LoopTask : public Task {
public:
    LoopTask() : Task(), functor_(nullptr) {}
    /*! \brief
     * Set the functor (work) to be done.
     *
     * Note: Takes ownership of passed functor
     *
     * This releases a barrier so that threads who have or will call run can
     * proceed.
     *
     * This function is not thread safe and is intended to be called only by
     * the main thread (thread running the containing basic task or thread 0
     * for the default schedule).
     *
     * \param[in] f Task functor
     */
    void setFunctor(ITaskFunctor* f) override {
        functor_.reset(f);
        functorBeginBarrier_.open();
    }
    /*! \brief
     * Return whether task is ready to run
     */
    bool isReady() const override {
        return functorBeginBarrier_.isOpen();
    }
    /*! \brief
     * Run the functor for the specified range
     *
     * This function is thread-safe and intended to be called by all
     * participating threads. It synchronizes threads and marks when they
     * complete.
     *
     * \param[in] range Range of function to run
     * \return whether all functors have been assigned for this task, which
     *         is always true for a LoopTask after running its single task.
     */
    bool run(Range<Ratio> range) override {
        functorBeginBarrier_.wait();
        functor_->run(range);
        functorEndBarrier_.markArrival();
        return true;
    }
    /*! \brief
     * Wait for all participating threads to finish
     *
     * Normally called by main thread but can be called by any thread to
     * wait for task to complete. Is thread-safe.
     */
    void wait() override {
        functorEndBarrier_.wait();
    }
private:
    /*! \brief
     * Reset this object for running a new functor. Nullifies any stored
     * functors and resets barriers. Intended only to be called by thread 0
     * in-between steps.
     */
    void init() override {
        functor_.reset(nullptr);
        functorBeginBarrier_.close();
        functorEndBarrier_.close(this->getNumThreads());
    }
    std::unique_ptr<ITaskFunctor> functor_;      //!< The function/loop to execute
    MOBarrier functorBeginBarrier_; //!< Many-to-one barrier to sync threads at beginning of loop
    OMBarrier functorEndBarrier_; //!< One-to-many barrier to sync threads at end of loop
};

class MultiLoopTask : public Task {
public:
    MultiLoopTask() : Task(), functor_(nullptr),
    functorCounter_(0) {}
    /*! \brief
     * Set the functor (work) to be done.
     *
     * Note: Takes ownership of passed functor
     *
     * Unlike LoopTask, an atomic counter is incremented when a new function is
     * assigned. More than a simple barrier is necessary to support multiple
     * functor assignments within a single step.
     *
     * This function is not thread safe and is intended to be called only by
     * the main thread (thread running the containing basic task). MultiLoops
     * do not occur with a default schedule.
     *
     * \param[in] f Task functor
     */
    void setFunctor(ITaskFunctor* f) override {
        functorEndBarrier_.close(this->getNumThreads());
        functor_.reset(f);
        functorCounter_++;
    }
    /*! \brief
     * Return whether task is ready to run
     */
    bool isReady() const override {
        int localThreadId = this->getThreadId(Thread::getId());
        return functorCounter_.load() != threadCounters_[localThreadId];
    }
    /*! \brief
     * Mark that no more functions will be assigned during the current step.
     *
     * While thread-safe, this is intended to only be called by the main thread
     * once the containing basic task completes.
     */
    void markFinished() {
        functorCounter_ = -1;
    }
    /*! \brief
     * Run the functor for the specified range
     *
     * This function is thread-safe and intended to be called by all
     * participating threads. It synchronizes threads and marks when they
     * complete. Unlike LoopTask, each thread has a unique counter used to
     * synchronize running multiple functors.
     *
     * \param[in] range Range of function to run
     * \return whether all functors have been assigned for this task. If
     *         not, thread should call run again.
     */
    bool run(Range<Ratio> range) override {
        int localThreadId = this->getThreadId(Thread::getId());
        wait_until_not(functorCounter_, threadCounters_[localThreadId]++);
        if (functorCounter_ == -1) {
            return true;
        }
        functor_->run(range);
        functorEndBarrier_.markArrival();
        return false;
    }
    /*! \brief
     * Wait for all participating threads to finish
     *
     * Normally called by main thread but can be called by any thread to
     * wait for task to complete. Is thread-safe.
     */
    void wait() override {
        functorEndBarrier_.wait();
    }
private:
    /*! \brief
     * Reset this object for running a new set of functors. Nullifies any
     * stored functors and resets counters and end barrier. Intended only to
     * be called by thread 0 in-between steps.
     */
    void init() override {
        functorCounter_ = 0;
        functor_.reset(nullptr);
        functorEndBarrier_.close(this->getNumThreads());
        threadCounters_.assign(this->getNumThreads(),0);
    }
    std::unique_ptr<ITaskFunctor> functor_;      //!< The function/loop to execute
    OMBarrier functorEndBarrier_; //!< One-to-many barrier to sync threads at end of loop
    std::atomic<int> functorCounter_;
    std::vector<int> threadCounters_;
};

class BasicTask : public Task {
public:
    BasicTask() : Task(), functor_(nullptr), containedMultiLoopTask_(nullptr) {}
    /*! \brief
     * Set the functor (work) to be done.
     *
     * Note: Takes ownership of passed functor
     *
     * This releases a barrier so that thread who has or will call run can
     * proceed.
     *
     * This function is not thread safe and is intended to be called only by
     * thread 0 (basic tasks cannot be nested in other tasks).
     *
     * \param[in] f Task functor
     */
    void setFunctor(ITaskFunctor* f) override {
        functor_.reset(f);
        functorBeginBarrier_.open();
    }
    /*! \brief
     * Return whether task is ready to run
     */
    bool isReady() const override {
        return functorBeginBarrier_.isOpen();
    }
    /*! \brief
     * Run the functor. The range argument is ignored.
     *
     * This function is thread-safe and intended to be called by the thread
     * assigned to this task. Thread waits until functor is available.
     *
     * \param[in] range Range to run - ignored.
     * \return whether all functors have been assigned for this task, which
     *         is always true for a BasicTask after running its single task.
     */
    bool run(Range<Ratio> range) override {
        functorBeginBarrier_.wait();
        functor_->run(range);
        functorEndBarrier_.markArrival();
        if (containedMultiLoopTask_ != nullptr) {
            containedMultiLoopTask_->markFinished();
        }
        return true;
    }
    /*! \brief
     * Set that the given multiloop task is contained within this task. This is
     * necessary so that this task can mark the multiloop task as finished once
     * it itself completes.
     *
     * \param[in] mlt MultiLoopTask contained in this task.
     */
    void setMultiLoop(MultiLoopTask *mlt) {
        assert(containedMultiLoopTask_ == nullptr);
        containedMultiLoopTask_ = mlt;
    }
    /*! \brief
     * Wait for all participating threads to finish
     *
     * Normally called by main thread but can be called by any thread to
     * wait for task to complete. Is thread-safe.
     */
    void wait() override {
        functorEndBarrier_.wait();
    }
private:
    /*! \brief
     * Reset this object for running a new functor. Nullifies any stored
     * functors and resets barriers. Intended only to be called by thread 0
     * in-between steps.
     */
    void init() override {
        functor_.reset(nullptr);
        functorBeginBarrier_.close();
        functorEndBarrier_.close(this->getNumThreads());
    }
    std::unique_ptr<ITaskFunctor> functor_;      //!< The function/loop to execute
    MOBarrier functorBeginBarrier_; //!< Many-to-one barrier to sync threads at beginning of loop
    OMBarrier functorEndBarrier_; //!< One-to-many barrier to sync threads at end of loop
    MultiLoopTask* containedMultiLoopTask_;
};

#endif // STS_TASK_H
