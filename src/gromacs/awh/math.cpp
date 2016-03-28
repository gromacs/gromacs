/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "math.h"

#include <assert.h>

#include <array>
#include <random>

#include "gromacs/math/utilities.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

int Math::getSampleFromDistribution(ArrayRef<const double> distr,
                                    gmx_int64_t            seed,
                                    gmx_int64_t            indexSeed0,
                                    gmx_int64_t            indexSeed1)
{
    gmx::ThreeFry2x64<0>               rng(seed, gmx::RandomDomain::AwhBiasing);
    gmx::UniformRealDistribution<real> uniformRealDistr;

    GMX_RELEASE_ASSERT(distr.size() > 0, "We need a non-zero length distribution to sample from");

    /* Generate the cumulative probability distribution function */
    std::vector<double> distrCumul(distr.size());

    distrCumul[0] = distr[0];

    for (size_t i = 1; i < distr.size(); i++)
    {
        distrCumul[i] = distrCumul[i - 1] + distr[i];
    }

    GMX_RELEASE_ASSERT(gmx_within_tol(distrCumul.back(), 1.0, 0.01), "Attempt to get sample from non-normalized/zero distribution");

    /* Use binary search to convert the real value to an integer in [0, ndistr - 1] distributed according to distr. */
    rng.restart(indexSeed0, indexSeed1);

    double value  = uniformRealDistr(rng);
    int    sample = std::upper_bound(distrCumul.begin(), distrCumul.end() - 1, value) - distrCumul.begin();

    return sample;
}

double Math::expSum(double a, double b)
{
    return (a > b ? a : b) + std::log1p(std::exp(-std::fabs(a - b)));
}

double Math::gaussianGeometryFactor(gmx::ArrayRef<const double> xArray)
{
    /* For convenience we give the geometry factor function a name: zeta(x) */
    constexpr size_t                    tableSize   = 5;
    std::array<const double, tableSize> xTabulated  =
    { {1e-5, 1e-4, 1e-3, 1e-2, 1e-1} };
    std::array<const double, tableSize> zetaTable1d =
    { { 0.166536811948, 0.16653116886, 0.166250075882, 0.162701098306, 0.129272430287 } };
    std::array<const double, tableSize> zetaTable2d =
    { { 2.31985974274, 1.86307292523, 1.38159772648, 0.897554759158, 0.405578211115 } };

    gmx::ArrayRef<const double> zetaTable;

    if (xArray.size() == 1)
    {
        zetaTable = zetaTable1d;
    }
    else if (xArray.size() == 2)
    {
        zetaTable = zetaTable2d;
    }
    else
    {
        /* TODO... but this is anyway a rough estimate and > 2 dimensions is not so popular. */
        zetaTable = zetaTable2d;
    }

    /* TODO. Really zeta is a function of an ndim-dimensional vector x and we shoudl have a ndim-dimensional lookup-table.
       Here we take the geometric average of the components of x which is ok if the x-components are not very different. */
    double xScalar = 1;
    for (auto &x : xArray)
    {
        xScalar *= x;
    }

    xScalar = std::pow(xScalar, 1.0/xArray.size());

    /* Look up zeta(x) */
    size_t xIndex = 0;
    while ((xIndex < xTabulated.size()) && (xScalar > xTabulated[xIndex]))
    {
        xIndex++;
    }

    double zEstimate;
    if (xIndex == xTabulated.size())
    {
        /* Take last value */
        zEstimate = zetaTable[xTabulated.size() - 1];
    }
    else if (xIndex == 0)
    {
        zEstimate = zetaTable[xIndex];
    }
    else
    {
        /* Interpolate */
        double x0 = xTabulated[xIndex - 1];
        double x1 = xTabulated[xIndex];
        double w  = (xScalar - x0)/(x1 - x0);
        zEstimate = w*zetaTable[xIndex - 1] + (1 - w)*zetaTable[xIndex];
    }

    return zEstimate;
}

} // namespace gmx
