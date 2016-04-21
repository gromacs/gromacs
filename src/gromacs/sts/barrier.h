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
template <typename A, typename T>
void wait_until(const A &a, T v) {
    while (a.load() != v);
}
/*! \brief
 * Wait until atomic variable a is not set to value v
 *
 * \param[in] a   atomic variable
 * \param[in] v   value
 * \returns       new value of a
 */
template <typename A, typename T>
T wait_until_not(const A &a, T v) {
    T v2;
    do {
        v2 = a.load();
    } while(v2 == v);
    return v2;
}

// TODO: Asserts to check for wrong or multiple uses of barriers.

/*! \brief
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
     * Reset barrier
     */
    void close() {isLocked.store(true);}
private:
    std::atomic_bool isLocked;
};

/*! \brief
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
    std::atomic_int numThreadsRemaining;
};

/*! \brief
 * A simple many-to-many (MM) barrier.
 *
 * This is a reusable barrier and so works inside loops.
 * It assumes a fixed set of exactly nt threads.
 */
class MMBarrier {
public:
    MMBarrier(int nt, std::string name = "") :id(name), nthreads(nt),
                                            numWaitingThreads(0),
                                            numReleasedThreads(0) {
        assert(nt > 0);
        if (!id.empty()) {
            barrierInstances_[id] = this;
        }
    }
    void enter() {
        wait_until(numReleasedThreads, 0);
        numWaitingThreads.fetch_add(1);
        wait_until(numWaitingThreads, nthreads);
        if (numReleasedThreads.fetch_add(1) == nthreads-1) {
            numWaitingThreads.store(0);
            numReleasedThreads.store(0);
        }
    }
    std::string getId() {
        return id;
    }
    /*! \brief
     * Returns MMBarrier instance for a given id or nullptr if not found
     *
     * \param[in] MMBarrier instance id
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

