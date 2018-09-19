/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMX_MATH_CONTAINERMEASURE_H
#define GMX_MATH_CONTAINERMEASURE_H

#include <cmath>

#include <algorithm>
#include <numeric>

#include "gromacs/utility/real.h"

namespace gmx
{

namespace containermeasure
{

template <typename T> typename T::value_type sum(const T &vecData)
{
    return std::accumulate(vecData.cbegin(), vecData.cend(), 0.);
}

//!\brief comparison operator where NaN is smaller than any number
template <typename T>
struct lessSmallNaN
{
    bool operator() ( T x, T y)
    {
        if (std::isnan(x))
        {
            return true;
        }
        if (std::isnan(y))
        {
            return false;
        }
        return x < y;
    };
};

//!\brief comparison operator where NaN is larger than any number
template <typename T>
struct lessBigNaN
{
    bool operator() ( T x, T y)
    {
        if (std::isnan(x))
        {
            return false;
        }
        if (std::isnan(y))
        {
            return true;
        }
        return x < y;
    };
};

//!\brief comparison operator where NaN is larger than any number
template <typename T>
struct sumNoNaN
{
    bool operator() ( T sum, const T element)
    {
        if (std::isnan(element))
        {
            return sum;
        }
        return sum+element;
    };
};


/*! \brief Sum, ignoring NaN values */
template <typename T> typename T::value_type sumIgnoreNaN(const T &vecData)
{
    return std::accumulate(vecData.cbegin(), vecData.cend(), 0., sumNoNaN<typename T::value_type>());
}

/*! \brief The minimum grid data value, excluding NaNs*/
template <typename T> typename T::value_type minIgnoreNaN(const T &vecData)
{
    return *std::min_element(vecData.cbegin(), vecData.cend(), lessBigNaN<typename T::value_type>());
}

/*! \brief The maximum grid data value, excluding NaNs. */
template <typename T> typename T::value_type maxIgnoreNaN(const T &vecData)
{
    return *std::max_element(vecData.cbegin(), vecData.cend(), lessSmallNaN<typename T::value_type>());
}

/*! \brief The minimum grid data value. */
template <typename T> typename T::value_type min(const T &vecData)
{
    return *std::min_element(vecData.cbegin(), vecData.cend(), lessBigNaN<typename T::value_type>());
}
/*! \brief The maximum grid data value. */
template <typename T> typename T::value_type max(const T &vecData)
{
    return *std::max_element(vecData.cbegin(), vecData.cend());
}

/*! \brief The mean grid data value. */
template <typename T> double meanIgnoreNaN(const T &vecData)
{
    const auto sum = sumIgnoreNaN(vecData);
    const auto num = std::count_if(vecData.cbegin(), vecData.cend(),  [](const typename T::value_type x){return std::isfinite(x); });
    if (num != 0)
    {
        return sum / static_cast<double>(num);
    }
    else
    {
        return std::numeric_limits<typename T::value_type>::quiet_NaN();
    }
}

/*! \brief The mean grid data value. */
template <typename T> double mean(const T &vecData)
{
    return sum(vecData) / static_cast<double>(vecData.size());
}

template <typename T> double sumOfSquareDeviationIgnoreNaN(const T &vecData)
{
    const auto data_mean = meanIgnoreNaN(vecData);
    return std::accumulate(
            vecData.cbegin(), vecData.cend(), 0.,
            [data_mean](const typename T::value_type &sum, typename T::value_type element) {
                if (std::isnan(element))
                {
                    return static_cast<double>(sum);
                }
                return sum + (element - data_mean)*(element - data_mean);
            });
}


template <typename T> double sumOfSquareDeviation(const T &vecData)
{
    double data_mean = mean(vecData);
    return std::accumulate(
            vecData.cbegin(), vecData.cend(), 0.,
            [data_mean](const typename T::value_type &a, typename T::value_type b) {
                return a + (b - data_mean)*(b - data_mean);
            });
}
template <typename T> double normSquared(const T &vecData)
{
    return std::accumulate(
            vecData.cbegin(), vecData.cend(), 0.,
            [](const typename T::value_type &a, const typename T::value_type &b) {
                return a + b*b;
            });
}

template <typename T> double norm(const T &vecData)
{
    return sqrt(normSquared(vecData));
}
/*! \brief The variance of data values.
 *
 * Two pass algorithm, where mean is calculated fist, due to large expected
 * amount of data. (Typical data size=1,000,000)
 */
template <typename T> double var(const T &vecData)
{
    return sumOfSquareDeviation(vecData) / vecData.size();
}

/*! \brief The variance of data values ignoring NaN values.
 *
 */
template <typename T> double varIgnoreNaN(const T &vecData)
{
    const auto sum = sumOfSquareDeviationIgnoreNaN(vecData);
    const auto num = std::count_if(vecData.cbegin(), vecData.cend(), [](const typename T::value_type x){return std::isfinite(x); });

    if (num != 0)
    {
        return sum / static_cast<double>(num);
    }
    return std::numeric_limits<typename T::value_type>::quiet_NaN();
}

/*! \brief The root mean square deviation of grid data values. */
template <typename T> double rms(const T &vecData)
{
    return sqrt(var(vecData));
}

/*! \brief The root mean square deviation of grid data values. */
template <typename T> double rmsIgnoreNaN(const T &vecData)
{
    return sqrt(varIgnoreNaN(vecData));
}

template <typename T>
double entropy(const T &vecData, typename T::value_type unitCellVolume)
{
    auto p    = vecData.cbegin();
    int  size = vecData.size();
    real sum  = 0;
    for (int i = 0; i < size; ++i)
    {
        if (p[i] > 0)
        {
            sum += p[i] * log(p[i]);
        }
    }
    return -unitCellVolume * sum;
}

/*!\brief Returns new container, minued - subtrahend. */
template <typename T> T subtract(const T &subtrahend, const T &minuend)
{
    T result;
    std::transform(std::begin(minuend), std::end(minuend), std::begin(subtrahend),
                   std::back_inserter(result),
                   [](float min, float sub) { return min - sub; });
    return result;
}

/*! \brief Returns new container, log(denominator/nominator). */
template <typename T> T logRatio(const T &nominator, const T &denominator)
{
    T result;
    std::transform(std::begin(nominator), std::end(nominator),
                   std::begin(denominator), std::back_inserter(result),
                   [](typename T::value_type nom, typename T::value_type den) {
                       return log(nom/ den);
                   });
    return result;
}

//! \brief Returns new container, -log(input). */
template <typename T> T negLog(const T &input)
{
    T result;
    std::transform(std::begin(input), std::end(input), std::back_inserter(result),
                   [](typename T::value_type value) { return -log(value); });
    return result;
}
}

} // namespace gmx
#endif
