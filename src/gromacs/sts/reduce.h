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

#ifndef STS_REDUCE_H
#define STS_REDUCE_H

#include <vector>

/*!
 * Reduction class and implementation for basic data types
 *
 * Use template specialization to implement reductions for other
 * types.
 *
 * These classes are meant to be passed into parallel_for loops.
 * Threads should call collect to contribute values. STS calls reduce
 * at the end of the parallel_for to combine all of the contributed
 * values. Afterwards, call getResult() to read the final result.
 *
 * This class is thread safe when used through the interface provided
 * in "sts.h" rather than when used directly. This interface gives
 * each thread a unique slot for contributing values, and reduce is
 * only called by the main thread after the parallel_for ends and not
 * by the user.
 *
 * Custom implementations should be careful to ensure thread safety.
 */
template<typename T>
class TaskReduction {
public:
    TaskReduction(T init, int numThreads) :result(init) {
        values.resize(numThreads, init);
    }
    void collect(T a, size_t pos) {
        values[pos] += a;
    }
    // TODO: Allow user to provide a custom reduction function
    void reduce() {
        for (const T &i : values) {
            result += i;
        }
    }
    T getResult() {
        return result;
    }
private:
    std::vector<T> values;
    T result;
};

#endif // STS_REDUCE_H
