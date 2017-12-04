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
#ifndef GMX_MATH_REALFIELDMEASURE_H
#define GMX_MATH_REALFIELDMEASURE_H

#include <vector>
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/gmxomp.h"


namespace gmx
{

namespace realgriddatameasure
{
RVec center_of_mass(const GridDataReal3D &realfield);
}

namespace datavectormeasure
{
template <typename T>
typename T::value_type sum(const T &vecData)
{
    return std::accumulate(vecData.cbegin(), vecData.cend(), 0.);
}

/*! \brief The minimum grid data value. */
template <typename T>
typename T::value_type min(const T &vecData)
{
    return *std::min_element(vecData.cbegin(), vecData.cend());
}
/*! \brief The maximum grid data value. */
template <typename T>
typename T::value_type max(const T &vecData)
{
    return *std::max_element(vecData.cbegin(), vecData.cend());
}
/*! \brief The mean grid data value. */
template <typename T>
double mean(const T &vecData)
{
    return sum(vecData) / static_cast<double>(vecData.size());
}

template <typename T>
double sumOfSquareDeviation(const T &vecData)
{
    real data_mean = mean(vecData);
    return std::accumulate(vecData.cbegin(), vecData.cend(), 0.,
                           [data_mean](const real &a, real b) {
                               return a + gmx::square(b - data_mean);
                           });
}

template <typename T>
double normSquared(const T &vecData)
{
    return std::accumulate(
            vecData.cbegin(), vecData.cend(), 0.,
            [](const typename T::value_type &a, const typename T::value_type &b) { return a + square(b); });
}

template <typename T>
double norm(const T &vecData) { return sqrt(normSquared(vecData)); }
/*! \brief The variance of data values.
 *
 * Two pass algorithm, where mean is calculated fist, due to large expected
 * amount of data. (Typical data size=1,000,000)
 */
template <typename T>
double var(const T &vecData)
{
    return sumOfSquareDeviation(vecData) / vecData.size();
}
/*! \brief The root mean square deviation of grid data values. */
template <typename T>
double rms(const T &vecData) { return sqrt(var(vecData)); }

template <typename T>
double entropy(const T &realfield, typename T::value_type unitCellVolume)
{
    auto  p    = realfield.cbegin();
    int   size = realfield.size();
    real  sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_get_max_threads())) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_get_max_threads()) + 1)
    for (int i = 0; i < size; ++i)
    {
        if (p[i] > 0)
        {
            sum += p[i] * log(p[i]);
        }
    }
    return -unitCellVolume * sum;
}
}

}      // namespace gmx
#endif
