/*
 * LogicallyErrorComparer.h
 *
 *  Created on: Oct 9, 2014
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef LOGICALLYERRORCOMPARER_H_
#define LOGICALLYERRORCOMPARER_H_

#include <cmath>
#include <cstdlib>
#include "gromacs/utility/real.h"

struct LogicallyEqualComparerBase
{
    LogicallyEqualComparerBase(double error_factor) : error_factor(error_factor) {}

    double error_factor;
};

template <bool weight_by_magnitude = true, bool ignore_sign = false>
struct LogicallyEqualComparer : public LogicallyEqualComparerBase
{
    LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

    template <class T>
    bool operator () (T a, T b) const
    {
        const T absA = std::abs(a);
        const T absB = std::abs(b);
        const T diff = std::abs(a-b);

        if (a == b) {
            return true;
        } else if ((absA + absB) < 1.0) {
            return diff < std::numeric_limits<real>::epsilon() * error_factor;
        } else {
            return diff < (absA + absB) * std::numeric_limits<real>::epsilon() * error_factor;
        }
    }
};

template <>
struct LogicallyEqualComparer<true,true> : public LogicallyEqualComparerBase
{
    LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

    template <class T>
    bool operator () (T a, T b) const
    {
        const T absA = std::abs(a);
        const T absB = std::abs(b);
        const T diff = absA - absB;

        if (a == b) {
            return true;
        } else if ((absA + absB) < 1.0) {
            return diff < std::numeric_limits<real>::epsilon() * error_factor;
        } else {
            return diff < (absA + absB) * std::numeric_limits<real>::epsilon() * error_factor;
        }
    }
};

template <>
struct LogicallyEqualComparer<false,false> : public LogicallyEqualComparerBase
{
    LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

    template <class T>
    bool operator () (T a, T b) const
    {
        const T absA = std::abs(a);
        const T absB = std::abs(b);
        const T diff = std::abs(a-b);

        if (a == b) {
            return true;
        } else {
            return diff < std::numeric_limits<real>::epsilon() * error_factor;
        }
    }
};

template <>
struct LogicallyEqualComparer<false,true> : public LogicallyEqualComparerBase
{
    LogicallyEqualComparer(double error_factor) : LogicallyEqualComparerBase(error_factor) {}

    template <class T>
    bool operator () (T a, T b) const
    {
        const T absA = std::abs(a);
        const T absB = std::abs(b);
        const T diff = absA - absB;

        if (a == b) {
            return true;
        } else {
            return diff < std::numeric_limits<real>::epsilon() * error_factor;
        }
    }
};

#endif /* LOGICALLYERRORCOMPARER_H_ */
