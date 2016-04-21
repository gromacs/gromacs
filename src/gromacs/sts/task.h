#ifndef STS_TASK_H
#define STS_TASK_H

#include <map>
#include <set>
#include <vector>

#include <chrono>
#include <iostream>

#include "barrier.h"
#include "range.h"
#include "thread.h"

using namespace std::chrono;
using sts_clock = steady_clock;
static const auto STS_MAX_TIME_POINT = time_point<sts_clock>::max();

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
template<typename F>
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

//! \internal Allows creating a loop functor without knowing the function type
template<typename F>
LoopTaskFunctor<F>* createLoopTaskFunctor(F f, Range<int64_t> r) {
    return new LoopTaskFunctor<F>(f,r);
}

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

//! \internal Allows creating a basic functor without knowing the function type
template<typename F>
BasicTaskFunctor<F>* createBasicTaskFunctor(F f) {
    return new BasicTaskFunctor<F>(f);
}

class Task;

struct TaskTimes {
public:
    time_point<sts_clock> waitStart;    // Time when work requested
    time_point<sts_clock> runStart;     // Time when work started
    time_point<sts_clock> runEnd;       // Time when work finished
    time_point<sts_clock> nextRunAvail; // Time when next run (subtask) was ready
    std::map< std::string, std::vector<time_point<sts_clock>> > auxTimes;
    void clear() {
        waitStart    = STS_MAX_TIME_POINT;
        runStart     = STS_MAX_TIME_POINT;
        runEnd       = STS_MAX_TIME_POINT;
        nextRunAvail = STS_MAX_TIME_POINT;
        auxTimes.clear();
    }
    TaskTimes() {
       clear();
    }
};

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
     * \param[in] tid      Id of thread assigned the subtask
     * \param[in] task     The task this is part of.
     * \param[in] range    Out of a possible range from 0 to 1, the section in
                           this part. Ignored for basic tasks.
     */
    SubTask(int tid, Task *task, Range<Ratio> range) :threadId_(tid),
    task_(task), range_(range), isDone_(false) {}
    /*! \brief
     * Run the subtask
     */
    bool run();
    const Task *getTask() const;
    bool isDone() const;
    void setDone(bool isDone);
    bool isReady() const;
    long getWaitStartTime() const {
        auto s = time_point_cast<microseconds>(timeData_.waitStart);
        return s.time_since_epoch().count();
    }
    long getRunStartTime() const {
        auto s = time_point_cast<microseconds>(timeData_.runStart);
        return s.time_since_epoch().count();
    }
    long getRunEndTime() const {
        auto s = time_point_cast<microseconds>(timeData_.runEnd);
        return s.time_since_epoch().count();
    }
    long getNextRunAvailTime() const {
        auto s = time_point_cast<microseconds>(timeData_.nextRunAvail);
        return s.time_since_epoch().count();
    }
    /*! \brief
     * Set time when the next subtask became available (lambda was set)
     *
     * \param[in] t time in microseconds since epoch
     */
    void setNextRunAvailTime(long t) {
       timeData_.nextRunAvail = time_point<sts_clock>(microseconds(t));
    }
    long getWaitDuration() const {
        return getRunStartTime() - getWaitStartTime();
    }
    long getRunDuration() const {
        return getRunEndTime() - getRunStartTime();
    }
    long getTotalDuration() const {
        return getRunEndTime() - getWaitStartTime();
    }
    void recordTime(std::string label) {
        timeData_.auxTimes[label].push_back(sts_clock::now());
    }
    std::vector<long> getAuxTimes(std::string label) const {
        std::vector<long> times;
        for (auto t : timeData_.auxTimes.at(label)) {
            auto s = time_point_cast<microseconds>(t);
            times.push_back(s.time_since_epoch().count());
        }
        return times;
    }
    void clearTimes() {
        timeData_.clear();
    }
    Range<Ratio> getRange() const {
        return range_;
    }
    void setRange(Range<Ratio> r) {
        range_ = r;
    }
    const int threadId_;
private:
    Task *task_;             /**< Reference to main task */
    Range<Ratio> range_;     /**< Range (out of [0,1]) of loop part */
    bool isDone_;
    TaskTimes timeData_;
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
    Task(std::string l) :reduction_(nullptr), label(l), numThreads_(0),
                         priority_(NORMAL), functorSetTime_(STS_MAX_TIME_POINT) {}
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
    //! \brief Get total number of subtasks for this task
    int getNumSubtasks() const {
        return subtasks_.size();
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
            st->clearTimes();
        }
        init();
    }
    //! \brief Get task label
    std::string getLabel() const {
        return label;
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
    const SubTask* getSubTask(size_t i) const {
        if (i < subtasks_.size()) {
            return subtasks_[i].get();
        }
        return nullptr;
    }
    std::vector<const SubTask*> getSubTasks() const {
        std::vector<const SubTask*> v;
        for (auto const& st : subtasks_) {
            v.push_back(st.get());
        }
        return v;
    }
    /*! \brief
     * Set ranges for the subtasks of a task
     *
     * \param[in] intervals vector of ratios marking start and end points for each subtask
     * Example: setSubTaskRanges("reduce", {0,{1,6},{3,6},{4,6},1}
     */
    void setSubTaskRanges(std::vector<Ratio> intervals) {
        assert(intervals.size() == subtasks_.size()+1);
        assert(intervals[0] == 0);
        assert(intervals.back() == 1);
        for (size_t i=0; i<subtasks_.size(); i++) {
            assert(intervals[i] <= intervals[i+1]);
            subtasks_[i]->setRange({intervals[i],intervals[i+1]});
        }
    }
    /*! \brief
     * Set the functor (work) to be done.
     *
     * In the base class, we simply log the time of setting the functor
     * (important for load balancing) and call the subclass implementation.
     *
     * Note: Implementatons assume ownership of the passed functor.
     *
     * \param[in] f Task functor
     */
    void setFunctor(ITaskFunctor* f) {
        setFunctorImpl(f);
        functorSetTime_ = sts_clock::now();
    }
    long getFunctorSetTime() const {
        auto s = time_point_cast<microseconds>(functorSetTime_);
        return s.time_since_epoch().count();
    }
    //! \brief Return whether task is ready to run
    virtual bool isReady() const = 0;
    //! \brief Run the functor in the specified range
    virtual bool run(Range<Ratio>, TaskTimes &td) = 0;
    //! \brief Wait for the functor to complete
    virtual void wait() = 0;
    virtual ~Task() {}
protected:
    void*    reduction_;
private:
    virtual void setFunctorImpl(ITaskFunctor*) = 0;
    std::string label;
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
    time_point<sts_clock> functorSetTime_;
};

/*! \brief
 * Class for loop tasks
 */
class LoopTask : public Task {
public:
    LoopTask(std::string label) : Task(label), functor_(nullptr) {}
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
    bool run(Range<Ratio> range, TaskTimes &td) override {
        td.waitStart = sts_clock::now();
        functorBeginBarrier_.wait();
        td.runStart = sts_clock::now();
        functor_->run(range);
        td.runEnd = sts_clock::now();
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
    void setFunctorImpl(ITaskFunctor* f) override {
        functor_.reset(f);
        functorBeginBarrier_.open();
    }
    /*! \brief
     * Reset this object for running a new functor. Nullifies any stored
     * functors and resets barriers. Intended only to be called by thread 0
     * in-between steps.
     */
    void init() override {
        functor_.reset(nullptr);
        functorBeginBarrier_.close();
        functorEndBarrier_.close(this->getNumSubtasks());
    }
    std::unique_ptr<ITaskFunctor> functor_;      //!< The function/loop to execute
    MOBarrier functorBeginBarrier_; //!< Many-to-one barrier to sync threads at beginning of loop
    OMBarrier functorEndBarrier_; //!< One-to-many barrier to sync threads at end of loop
};

class MultiLoopTask : public Task {
public:
    MultiLoopTask(std::string label) : Task(label), functor_(nullptr),
    functorCounter_(0) {}
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
    bool run(Range<Ratio> range, TaskTimes &td) override {
        int localThreadId = this->getThreadId(Thread::getId());
        td.waitStart = sts_clock::now();
        wait_until_not(functorCounter_, threadCounters_[localThreadId]++);
        if (functorCounter_ == -1) {
            return true;
        }
        td.runStart = sts_clock::now();
        functor_->run(range);
        td.runEnd = sts_clock::now();
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
    void setFunctorImpl(ITaskFunctor* f) override {
        functorEndBarrier_.close(this->getNumThreads());
        functor_.reset(f);
        functorCounter_++;
    }
    /*! \brief
     * Reset this object for running a new set of functors. Nullifies any
     * stored functors and resets counters and end barrier. Intended only to
     * be called by thread 0 in-between steps.
     */
    void init() override {
        functorCounter_ = 0;
        functor_.reset(nullptr);
        functorEndBarrier_.close(this->getNumSubtasks());
        threadCounters_.assign(this->getNumThreads(),0);
    }
    std::unique_ptr<ITaskFunctor> functor_;      //!< The function/loop to execute
    OMBarrier functorEndBarrier_; //!< One-to-many barrier to sync threads at end of loop
    std::atomic<int> functorCounter_;
    std::vector<int> threadCounters_;
};

class BasicTask : public Task {
public:
    BasicTask(std::string label) : Task(label), functor_(nullptr), containedMultiLoopTask_(nullptr) {}
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
    bool run(Range<Ratio> range, TaskTimes &td) override {
        td.waitStart = sts_clock::now();
        functorBeginBarrier_.wait();
        td.runStart = sts_clock::now();
        functor_->run(range);
        td.runEnd = sts_clock::now();
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
    void setFunctorImpl(ITaskFunctor* f) override {
        functor_.reset(f);
        functorBeginBarrier_.open();
    }
    /*! \brief
     * Reset this object for running a new functor. Nullifies any stored
     * functors and resets barriers. Intended only to be called by thread 0
     * in-between steps.
     */
    void init() override {
        functor_.reset(nullptr);
        functorBeginBarrier_.close();
        functorEndBarrier_.close(this->getNumSubtasks());
    }
    std::unique_ptr<ITaskFunctor> functor_;      //!< The function/loop to execute
    MOBarrier functorBeginBarrier_; //!< Many-to-one barrier to sync threads at beginning of loop
    OMBarrier functorEndBarrier_; //!< One-to-many barrier to sync threads at end of loop
    MultiLoopTask* containedMultiLoopTask_;
};

#endif // STS_TASK_H
