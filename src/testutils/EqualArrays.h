/*
 * EqualArrays.h
 *
 *  Created on: Feb 20, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef EQUALARRAYS_H_
#define EQUALARRAYS_H_

#include "gtest/gtest.h"
#include <cmath>
#include <vector>

template <class T>
struct DefaultTolerance
{
	static const T value = 0;
};

//! Threshold for equality check of two floating point numbers.
template <>
struct DefaultTolerance<float>
{
	static constexpr float value = 1.0 / (1 << 22);
};

//! Threshold for equality check of two floating point numbers.
template <>
struct DefaultTolerance<double>
{
	static constexpr double value = 1.0 / (1ul << 51);
};

//! Check equality of two arrays.
template <class T>
::testing::AssertionResult EqualArrays(const T* const expected,
    const T* const actual, unsigned long length, float tolerance = DefaultTolerance<T>::value)
{
    ::testing::AssertionResult result = ::testing::AssertionFailure();
    int errorsFound = 0;
    const char* separator = " ";
    for (unsigned long index = 0; index < length; index++)
    {
        if (std::abs(expected[index] - actual[index]) > tolerance)
        {
            if (errorsFound == 0)
            {
                result << "Differences found:";
            }
            if (errorsFound < 3)
            {
                result << separator
                        << expected[index] << " != " << actual[index]
                        << " @ " << index;
                separator = ", ";
            }
            errorsFound++;
        }
    }
    if (errorsFound > 0)
    {
        result << separator << errorsFound << " differences in total";
        result << separator << "tolerance = " << tolerance;
        return result;
    }
    return ::testing::AssertionSuccess();
}

//! Check equality of two std::vectors.
template <class T>
::testing::AssertionResult EqualArrays(std::vector<T> const& expected,
	std::vector<T> const& actual, T tolerance = DefaultTolerance<T>::value)
{
	if (actual.size() != expected.size()) {
	    ::testing::AssertionResult result = ::testing::AssertionFailure();
	    result << "Difference in length: " << expected.size() << " != " << actual.size();
	    return result;
	}
	return EqualArrays(&expected[0], &actual[0], actual.size(), tolerance);
}

#endif /* EQUALARRAYS_H_ */
