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

#include "comparefields.h"
#include "fouriertransform.h"
#include "gromacs/math/griddata/griddata.h"
#include "gridinterpolator.h"
#include "realfieldmeasure.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include <map>
#include <iterator>
#include <set>
#include "gromacs/utility/exceptions.h"
#include <numeric>

namespace gmx
{

ComapreGridData::ComapreGridData(const GridDataReal3D &reference, const GridDataReal3D &other)
    : reference_ {reference}, other_ {
    other
}
{};

real ComapreGridData::correlate_(const std::vector<real> &a,
                                 const std::vector<real> &b) const
{

    auto              aMean = std::accumulate(a.cbegin(), a.cend(), 0.) / a.size();
    auto              aSSE  = std::accumulate(a.cbegin(), a.cend(), 0., [ aMean ](const real &a, real b) { return a + gmx::square(b - aMean); });
    auto              bMean = std::accumulate(b.cbegin(), b.cend(), 0.) / b.size();
    auto              bSSE  = std::accumulate(b.cbegin(), b.cend(), 0., [ bMean ](const real &a, real b) { return a + gmx::square(b - bMean); });
    std::vector<real> mulArray;
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(mulArray), [aMean, bMean](real aValue, real bValue) { return (aValue - aMean) * (bValue - bMean); });
    return std::accumulate(mulArray.begin(), mulArray.end(), 0.) / sqrt(aSSE*bSSE);
}

real ComapreGridData::correlate(real threshold) const
{
    std::vector<real> referenceAboveThreshold;
    std::vector<real> otherWhereReferenceAboveThreshold;

    auto              otherDatum = other_.begin();
    for (auto &referenceDatum : reference_)
    {
        if (referenceDatum > threshold)
        {
            referenceAboveThreshold.push_back(referenceDatum);
            otherWhereReferenceAboveThreshold.push_back(*otherDatum);
        }
        ++otherDatum;
    }

    return correlate_(referenceAboveThreshold, otherWhereReferenceAboveThreshold);
};

real ComapreGridData::gridSumAtCoordiantes(const std::vector<RVec> &coordinates)
{
    real        sum    = 0;
    for (const auto &r : coordinates)
    {
        sum += getLinearInterpolationAt(reference_, {{r[XX], r[YY], r[ZZ]}});
    }
    return sum;
};

real ComapreGridData::getRelativeKLCrossTermSameGrid(const std::vector<real> &other_reference) const
{
    // for numerical stability use a reference density
    if (reference_.size() != other_.size())
    {
        GMX_THROW(APIError(
                          "KL-Divergence calculation requires euqally sized input vectors."));
    }
    auto p           = reference_.begin();
    auto q           = other_.begin();
    auto q_reference = other_reference.begin();
    int  size        = other_.size();
    real sum         = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0) && (q_reference[i] > 0))
        {
            sum += p[i] * log(q[i] / (q_reference[i]));
        }
    }
    return -1 * reference_.getGrid().unitCell().volume() *sum;
}

real ComapreGridData::getKLSameGrid() const
{
    if (reference_.size() != other_.size())
    {
        GMX_THROW(APIError("KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p    = reference_.begin();
    auto q    = other_.begin();
    int  size = other_.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        real pp = p[i];
        if (pp > 0)
        {
            real r = q[i] / pp;
            if (r > 0)
            {
                sum += pp * log(r);
            }
        }

    }

    return -1 * reference_.getGrid().unitCell().volume() * sum;
};

real ComapreGridData::getKLCrossTermSameGrid() const
{
    if (reference_.size() != other_.size())
    {
        GMX_THROW(APIError(
                          "KL-CrossTerm calculation requires euqally sized input vectors."));
    }
    auto p    = reference_.begin();
    auto q    = other_.begin();
    int  size = other_.size();
    real sum  = 0;
#pragma omp parallel for num_threads(std::max( \
    1, gmx_omp_nthreads_get(emntDefault))) reduction(+ : sum) schedule( \
    static, size / std::max(1, gmx_omp_nthreads_get(emntDefault)) + 1)
    for (int i = 0; i < size; ++i)
    {
        if ((p[i] > 0) && (q[i] > 0))
        {
            sum += p[i] * log(q[i]);
        }
    }
    return -1 * reference_.getGrid().unitCell().volume() *sum;
};


}
