/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
#include <cmath>
#include <random>

#include "gromacs/math/utilities.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"


static int binary_search(double value, const double *array, int arraysize)
{
    int  imin  = 0;             /* min index */
    int  imax  = arraysize - 1; /* max index */
    int  imid;                  /* mid-point reference index */
    bool found = FALSE;
    do
    {
        /* start with (lower) middle element */
        imid = (imin + imax)/2;

        /* search in upper then lower array half */
        if (value >= array[imid])
        {
            imin = imid + 1;
        }
        else if (imid > 0 && value <  array[imid - 1])
        {
            imax = imid - 1;
        }
        else
        {
            found = TRUE;
        }

        /* Only one element left. We assume here that the value is in the array range. */
        if (imin == imax)
        {
            found = TRUE;
            imid  = imax;
        }
    }
    while (!found);

    return imid;
}

int get_sample_from_distribution(const double *distr, int ndistr, gmx_int64_t step, gmx_int64_t seed, int iseed)
{
    int                                sample;
    double                             value;
    double                            *distr_cumul;
    gmx::ThreeFry2x64<0>               rng(seed, gmx::RandomDomain::AwhBiasing);
    gmx::UniformRealDistribution<real> uniformRealDistr;

    GMX_RELEASE_ASSERT(ndistr > 0, "Attempt to get sample with from zero length distribution");

    /* The cumulative probability distribution function */
    snew(distr_cumul, ndistr);

    distr_cumul[0] = distr[0];

    for (int i = 1; i < ndistr; i++)
    {
        distr_cumul[i] = distr_cumul[i - 1] + distr[i];
    }

    GMX_RELEASE_ASSERT(distr_cumul[ndistr - 1] > 0.99, "Attempt to get sample from non-normalized/zero distribution");

    /* Use binary search to convert the real value to an integer in [0, ndistr - 1] distributed according to distr. */
    rng.restart(step, iseed);

    value  = uniformRealDistr(rng);
    sample = binary_search(value, distr_cumul, ndistr);

    sfree(distr_cumul);

    return sample;
}

double expsum(double a, double b)
{
    return (a > b ? a : b) + std::log1p(std::exp(-std::fabs(a - b)));
}

double gaussian_geometry_factor(const double *xarray, int ndim)
{
    /* For convenience we give the geometry factor function a name: zeta(x) */
    const double  x_tabulated[] =
    {1e-5, 1e-4, 1e-3, 1e-2, 1e-1};
    const double  zeta_table_1d[] =
    {
        0.166536811948, 0.16653116886, 0.166250075882,
        0.162701098306, 0.129272430287
    };
    const double  zeta_table_2d[] =
    {
        2.31985974274, 1.86307292523, 1.38159772648,
        0.897554759158, 0.405578211115
    };

    const double       *zeta_table;
    const int           npts_tabulated = asize(x_tabulated);
    double              x_scalar, z_estimate;

    if (ndim == 1)
    {
        zeta_table = zeta_table_1d;
    }
    else if (ndim == 2)
    {
        zeta_table = zeta_table_2d;
    }
    else
    {
        /* TODO... but this is anyway a rough estimate and > 2 dimensions is not so popular. */
        zeta_table = zeta_table_2d;
    }

    /* TODO. Really zeta is a function of an ndim-dimensional vector x and we shoudl have a ndim-dimensional lookup-table.
       Here we take the geometric average of the components of x which is ok if the x-components are not very different. */
    x_scalar = 1.;
    for (int d = 0; d < ndim; d++)
    {
        x_scalar *= xarray[d];
    }

    x_scalar = std::pow(x_scalar, 1./ndim);

    /* Look up zeta(x) */
    int xindex = 0;
    while ((xindex < npts_tabulated) && (x_scalar > x_tabulated[xindex]))
    {
        xindex++;
    }

    if (xindex == npts_tabulated)
    {
        xindex = npts_tabulated - 1;

        /* Take last value */
        z_estimate = zeta_table[xindex];
    }
    else if (xindex == 0)
    {
        z_estimate = zeta_table[xindex];
    }
    else
    {
        double x0, x1, w;

        /* Interpolate */
        x0          = x_tabulated[xindex - 1];
        x1          = x_tabulated[xindex];
        w           = (x_scalar - x0)/(x1 - x0);
        z_estimate  = w*zeta_table[xindex - 1] + (1 - w)*zeta_table[xindex];
    }
    return z_estimate;
}

int ceil_log2(double x)
{
    int k     = 0;
    int pow2k = 1;

    while (pow2k < x)
    {
        pow2k <<= 1;
        k++;
    }

    return k;
}
