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
    TaskReduction(std::string taskName, T init, int numThreads) :result(init) {
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
