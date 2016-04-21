#ifndef STS_BARRIER_H
#define STS_BARRIER_H

#include <cassert>

#include <map>
#include <string>

#include <atomic>

/*
 * Functions for spin waiting
 */
// TODO: Add gmx_pause.
// TODO: Consider wait_until_not variant that doesn't return a value.
/*! \brief
 * Wait until atomic variable a is set to value v
 *
 * \param[in] a   atomic variable
 * \param[in] v   value
 */
template <typename T>
void wait_until(const std::atomic<T> &a, T v) {
    while (a.load() != v);
}
//! Template overload to support nullptr
template <typename T>
void wait_until(const std::atomic<T> &a, std::nullptr_t) {
    while (a.load() != nullptr);
}
/*! \brief
 * Wait until atomic variable a is not set to value v
 *
 * \param[in] a   atomic variable
 * \param[in] v   value
 * \returns       new value of a
 */
template <typename T>
T wait_until_not(const std::atomic<T> &a, T v) {
    T v2;
    do {
        v2 = a.load();
    } while(v2 == v);
    return v2;
}
//! Template overload to support nullptr
template <typename T>
T wait_until_not(const std::atomic<T> &a, std::nullptr_t) {
    T v2;
    do {
        v2 = a.load();
    } while(v2 == nullptr);
    return v2;
}

// TODO: Asserts to check for wrong or multiple uses of barriers.

/*! \internal \brief
 * A simple many-to-one (MO) barrier.
 */
class MOBarrier {
public:
    MOBarrier() :isLocked(true) {}
    /*! \brief
     * Wait on barrier. Should be called by "M" threads
     */
    void wait() {
        wait_until(isLocked, false);
    }
    /*! \brief
     * Open barrier. Should be called by "O" thread.
     */
    void open() {isLocked.store(false);}
    /*! \brief
     * Check if barrier is open
     */
    bool isOpen() const {return !isLocked.load();}
    /*! \brief
     * Reset barrier
     */
    void close() {isLocked.store(true);}
private:
    std::atomic<bool> isLocked;
};

/*! \internal \brief
 * A simple one-to-many (OM) barrier.
 */
class OMBarrier {
public:
    OMBarrier() :numThreadsRemaining(0) {}
    /*! \brief
     * Register with the barrier. Should be called by "M" threads.
     */
    void markArrival() {
        numThreadsRemaining--;
    }
    /*! \brief
     * Wait on barrier. Should be called by "O" thread.
     */
    void wait() {
        wait_until(numThreadsRemaining, 0);
    }
    /*! \brief
     * Reset barrier
     */
    void close(int nthreads) {
        numThreadsRemaining.store(nthreads);
    }

private:
    std::atomic<int> numThreadsRemaining;
};

/*! \internal \brief
 * A simple many-to-many (MM) barrier.
 *
 * This is a reusable barrier and so works inside loops.
 * It assumes a fixed set of exactly nt threads.
 */
class MMBarrier {
public:
    /*! \brief
     * Constructs a new many-to-many barrier.
     *
     * \param[in] nt    Number of threads
     * \param[in] name  Barrier name
     */
    MMBarrier(int nt, std::string name = "") :id(name), nthreads(nt),
                                            numWaitingThreads(0),
                                            numReleasedThreads(0) {
        assert(nt > 0);
        if (!id.empty()) {
            barrierInstances_[id] = this;
        }
    }
    //! \brief Enter barrier
    void enter() {
        wait_until(numReleasedThreads, 0);
        numWaitingThreads.fetch_add(1);
        wait_until(numWaitingThreads, nthreads);
        if (numReleasedThreads.fetch_add(1) == nthreads-1) {
            numWaitingThreads.store(0);
            numReleasedThreads.store(0);
        }
    }
    /*! \brief
     * Get barrier id (name)
     *
     * \returns barrier id (name)
     */
    std::string getId() {
        return id;
    }
    /*! \brief
     * Returns MMBarrier instance for a given id or nullptr if not found
     *
     * \param[in] id  MMBarrier instance id
     * \returns MMBarrier instance
     */
    static MMBarrier *getInstance(std::string id) {
        auto entry = barrierInstances_.find(id);
        if (entry == barrierInstances_.end()) {
            return nullptr;
        }
        else {
            return entry->second;
        }
    }
private:
    std::string id;
    const int nthreads;
    std::atomic<int>  numWaitingThreads;
    std::atomic<int>  numReleasedThreads;
    static std::map<std::string, MMBarrier *> barrierInstances_;
};

#endif // STS_BARRIER_H

