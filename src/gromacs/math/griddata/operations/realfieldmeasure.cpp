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

#include "realfieldmeasure.h"
#include "gromacs/utility/gmxomp.h"

namespace gmx
{

DataVectorMeasure::DataVectorMeasure(const std::vector<float> &vecData_)
    : vecData_ {vecData_}
{};

float DataVectorMeasure::sum() const
{
    return std::accumulate(vecData_.cbegin(), vecData_.cend(), 0.);
}

/*! \brief The minimum grid data value. */
float DataVectorMeasure::min() const
{
    return *std::min_element(vecData_.cbegin(), vecData_.cend());
}
/*! \brief The maximum grid data value. */
float DataVectorMeasure::max() const
{
    return *std::max_element(vecData_.cbegin(), vecData_.cend());
}
/*! \brief The mean grid data value. */
float DataVectorMeasure::mean() const
{
    return sum() / static_cast<float>(vecData_.size());
}
/*! \brief
 * The size of the griddata in bytes.
 */
size_t DataVectorMeasure::data_size() const
{
    return vecData_.size();
}

float DataVectorMeasure::sumOfSquareDeviation() const
{
    float data_mean = mean();
    return std::accumulate(vecData_.cbegin(), vecData_.cend(), 0.,
                           [data_mean](const float &a, float b) {
                               return a + gmx::square(b - data_mean);
                           });
}

float DataVectorMeasure::normSquared() const
{
    return std::accumulate(
            vecData_.cbegin(), vecData_.cend(), 0.,
            [](const float &a, const float &b) { return a + square(b); });
}

float DataVectorMeasure::norm() const { return sqrt(normSquared()); }
/*! \brief The variance of data values.
 *
 * Two pass algorithm, where mean is calculated fist, due to large expected
 * amount of data. (Typical data size=1,000,000)
 */
float DataVectorMeasure::var() const
{
    return sumOfSquareDeviation() / vecData_.size();
}
/*! \brief The root mean square deviation of grid data values. */
float DataVectorMeasure::rms() const { return sqrt(var()); }

RealGridDataMeasure::RealGridDataMeasure(const GridDataReal3D &realfield)
    : realfield_ {realfield}
{};

RVec RealGridDataMeasure::center_of_mass() const
{
    RVec com = {0, 0, 0};
    for (size_t i = 0; i < realfield_.size(); i++)
    {
        auto multiIndex     = realfield_.getGrid().lattice().vectoriseLinearIndex(i);
        auto gridCoordinate = realfield_.getGrid().multiIndexToCoordinate(multiIndex);
        for (int dimension = XX; dimension <= ZZ; ++dimension)
        {
            com[dimension] += realfield_[i]*gridCoordinate[dimension];
        }
    }
    svmul(1. / (DataVectorMeasure(realfield_).sum()), com, com);
    return com;
}


float RealGridDataMeasure::entropy() const
{
    auto  p    = realfield_.cbegin();
    int   size = realfield_.size();
    float sum  = 0;
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
    return -1 * realfield_.getGrid().unitCell().volume() * sum;
};

} // namespace gmx
