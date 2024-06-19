/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "calculate_spline_moduli.h"

#include <cmath>

#include <algorithm>

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

static std::vector<real> make_dft_mod(gmx::ArrayRef<const double> data, int splineOrder, int ndata)
{
    std::vector<real> mod(ndata);

    for (int i = 0; i < ndata; i++)
    {
        /* We use double precision, since this is only called once per grid.
         * But for single precision bsp_mod, single precision also seems
         * to give full accuracy.
         */
        double sc = 0;
        double ss = 0;
        for (int j = 0; j < splineOrder; j++)
        {
            double arg = (2.0 * M_PI * i * (j + 1)) / ndata;
            sc += data[j] * std::cos(arg);
            ss += data[j] * std::sin(arg);
        }
        mod[i] = sc * sc + ss * ss;
    }
    if (splineOrder % 2 == 0 && ndata % 2 == 0)
    {
        /* Note that pme_order = splineOrder + 1 */
        GMX_RELEASE_ASSERT(mod[ndata / 2] < GMX_DOUBLE_EPS,
                           "With even spline order and even grid size (ndata), dft_mod[ndata/2] "
                           "should first come out as zero");
        /* This factor causes a division by zero. But since this occurs in
         * the tail of the distribution, the term with this factor can
         * be ignored (see Essmann et al. JCP 103, 8577).
         * Using the average of the neighbors probably originates from
         * Tom Darden's original PME code. It seems to give slighlty better
         * accuracy than using a large value.
         */
        mod[ndata / 2] = (mod[ndata / 2 - 1] + mod[ndata / 2 + 1]) * 0.5;
    }

    return mod;
}

std::array<std::vector<real>, 3> make_bspline_moduli(int nx, int ny, int nz, int pme_order)
{
    /* In GROMACS we, confusingly, defined pme-order as the order
     * of the cardinal B-spline + 1. This probably happened because
     * the smooth PME paper only talks about "n" which is the number
     * of points we spread to and that was chosen to be pme-order.
     */
    const int splineOrder = pme_order - 1;

    /* We use double precision, since this is only called once per grid.
     * But for single precision bsp_mod, single precision also seems
     * to give full accuracy.
     */
    std::vector<double> data(splineOrder);

    data[0] = 1;
    for (int k = 1; k < splineOrder; k++)
    {
        data[k] = 0;
    }

    for (int k = 2; k <= splineOrder; k++)
    {
        double div = 1.0 / k;
        for (int m = k - 1; m > 0; m--)
        {
            data[m] = div * ((k - m) * data[m - 1] + (m + 1) * data[m]);
        }
        data[0] = div * data[0];
    }

    std::array<std::vector<real>, 3> bsp_mod;

    bsp_mod[XX] = make_dft_mod(data, splineOrder, nx);
    bsp_mod[YY] = make_dft_mod(data, splineOrder, ny);
    bsp_mod[ZZ] = make_dft_mod(data, splineOrder, nz);

    return bsp_mod;
}

/* Return the P3M optimal influence function */
static double do_p3m_influence(double z, int order)
{
    double z2, z4;

    z2 = z * z;
    z4 = z2 * z2;

    /* The formula and most constants can be found in:
     * Ballenegger et al., JCTC 8, 936 (2012)
     */
    switch (order)
    {
        case 2: return 1.0 - 2.0 * z2 / 3.0;
        case 3: return 1.0 - z2 + 2.0 * z4 / 15.0;
        case 4: return 1.0 - 4.0 * z2 / 3.0 + 2.0 * z4 / 5.0 + 4.0 * z2 * z4 / 315.0;
        case 5:
            return 1.0 - 5.0 * z2 / 3.0 + 7.0 * z4 / 9.0 - 17.0 * z2 * z4 / 189.0 + 2.0 * z4 * z4 / 2835.0;
        case 6:
            return 1.0 - 2.0 * z2 + 19.0 * z4 / 15.0 - 256.0 * z2 * z4 / 945.0
                   + 62.0 * z4 * z4 / 4725.0 + 4.0 * z2 * z4 * z4 / 155925.0;
        case 7:
            return 1.0 - 7.0 * z2 / 3.0 + 28.0 * z4 / 15.0 - 16.0 * z2 * z4 / 27.0 + 26.0 * z4 * z4 / 405.0
                   - 2.0 * z2 * z4 * z4 / 1485.0 + 4.0 * z4 * z4 * z4 / 6081075.0;
        case 8:
            return 1.0 - 8.0 * z2 / 3.0 + 116.0 * z4 / 45.0 - 344.0 * z2 * z4 / 315.0
                   + 914.0 * z4 * z4 / 4725.0 - 248.0 * z4 * z4 * z2 / 22275.0
                   + 21844.0 * z4 * z4 * z4 / 212837625.0 - 8.0 * z4 * z4 * z4 * z2 / 638512875.0;
    }

    return 0.0;
}

/* Calculate the P3M B-spline moduli for one dimension */
static std::vector<real> make_p3m_bspline_moduli_dim(int n, int order)
{
    double zarg, zai, sinzai, infl;
    int    maxk, i;

    if (order > 8)
    {
        GMX_THROW(gmx::InconsistentInputError("The current P3M code only supports orders up to 8"));
    }

    std::vector<real> bsp_mod(n);

    zarg = M_PI / n;

    maxk = (n + 1) / 2;

    for (i = -maxk; i < 0; i++)
    {
        zai            = zarg * i;
        sinzai         = std::sin(zai);
        infl           = do_p3m_influence(sinzai, order);
        bsp_mod[n + i] = infl * infl * std::pow(sinzai / zai, -2.0 * order);
    }
    bsp_mod[0] = 1.0;
    for (i = 1; i < maxk; i++)
    {
        zai        = zarg * i;
        sinzai     = std::sin(zai);
        infl       = do_p3m_influence(sinzai, order);
        bsp_mod[i] = infl * infl * std::pow(sinzai / zai, -2.0 * order);
    }

    return bsp_mod;
}

/* Calculate the P3M B-spline moduli */
std::array<std::vector<real>, 3> make_p3m_bspline_moduli(int nx, int ny, int nz, int order)
{
    std::array<std::vector<real>, 3> bsp_mod;

    bsp_mod[XX] = make_p3m_bspline_moduli_dim(nx, order);
    bsp_mod[YY] = make_p3m_bspline_moduli_dim(ny, order);
    bsp_mod[ZZ] = make_p3m_bspline_moduli_dim(nz, order);

    return bsp_mod;
}
