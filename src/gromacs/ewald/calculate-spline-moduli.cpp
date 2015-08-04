/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "calculate-spline-moduli.h"

#include <math.h>

#include <algorithm>

#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"

static void make_dft_mod(real *mod, real *data, int ndata)
{
    int  i, j;
    real sc, ss, arg;

    for (i = 0; i < ndata; i++)
    {
        sc = ss = 0;
        for (j = 0; j < ndata; j++)
        {
            arg = (2.0*M_PI*i*j)/ndata;
            sc += data[j]*cos(arg);
            ss += data[j]*sin(arg);
        }
        mod[i] = sc*sc+ss*ss;
    }
    for (i = 0; i < ndata; i++)
    {
        if (mod[i] < 1e-7)
        {
            mod[i] = (mod[i-1]+mod[i+1])*0.5;
        }
    }
}

void make_bspline_moduli(splinevec bsp_mod,
                         int nx, int ny, int nz, int order)
{
    int   nmax = std::max(nx, std::max(ny, nz));
    real *data, *ddata, *bsp_data;
    int   i, k, l;
    real  div;

    snew(data, order);
    snew(ddata, order);
    snew(bsp_data, nmax);

    data[order-1] = 0;
    data[1]       = 0;
    data[0]       = 1;

    for (k = 3; k < order; k++)
    {
        div       = 1.0/(k-1.0);
        data[k-1] = 0;
        for (l = 1; l < (k-1); l++)
        {
            data[k-l-1] = div*(l*data[k-l-2]+(k-l)*data[k-l-1]);
        }
        data[0] = div*data[0];
    }
    /* differentiate */
    ddata[0] = -data[0];
    for (k = 1; k < order; k++)
    {
        ddata[k] = data[k-1]-data[k];
    }
    div           = 1.0/(order-1);
    data[order-1] = 0;
    for (l = 1; l < (order-1); l++)
    {
        data[order-l-1] = div*(l*data[order-l-2]+(order-l)*data[order-l-1]);
    }
    data[0] = div*data[0];

    for (i = 0; i < nmax; i++)
    {
        bsp_data[i] = 0;
    }
    for (i = 1; i <= order; i++)
    {
        bsp_data[i] = data[i-1];
    }

    make_dft_mod(bsp_mod[XX], bsp_data, nx);
    make_dft_mod(bsp_mod[YY], bsp_data, ny);
    make_dft_mod(bsp_mod[ZZ], bsp_data, nz);

    sfree(data);
    sfree(ddata);
    sfree(bsp_data);
}

/* Return the P3M optimal influence function */
static double do_p3m_influence(double z, int order)
{
    double z2, z4;

    z2 = z*z;
    z4 = z2*z2;

    /* The formula and most constants can be found in:
     * Ballenegger et al., JCTC 8, 936 (2012)
     */
    switch (order)
    {
        case 2:
            return 1.0 - 2.0*z2/3.0;
            break;
        case 3:
            return 1.0 - z2 + 2.0*z4/15.0;
            break;
        case 4:
            return 1.0 - 4.0*z2/3.0 + 2.0*z4/5.0 + 4.0*z2*z4/315.0;
            break;
        case 5:
            return 1.0 - 5.0*z2/3.0 + 7.0*z4/9.0 - 17.0*z2*z4/189.0 + 2.0*z4*z4/2835.0;
            break;
        case 6:
            return 1.0 - 2.0*z2 + 19.0*z4/15.0 - 256.0*z2*z4/945.0 + 62.0*z4*z4/4725.0 + 4.0*z2*z4*z4/155925.0;
            break;
        case 7:
            return 1.0 - 7.0*z2/3.0 + 28.0*z4/15.0 - 16.0*z2*z4/27.0 + 26.0*z4*z4/405.0 - 2.0*z2*z4*z4/1485.0 + 4.0*z4*z4*z4/6081075.0;
        case 8:
            return 1.0 - 8.0*z2/3.0 + 116.0*z4/45.0 - 344.0*z2*z4/315.0 + 914.0*z4*z4/4725.0 - 248.0*z4*z4*z2/22275.0 + 21844.0*z4*z4*z4/212837625.0 - 8.0*z4*z4*z4*z2/638512875.0;
            break;
    }

    return 0.0;
}

/* Calculate the P3M B-spline moduli for one dimension */
static void make_p3m_bspline_moduli_dim(real *bsp_mod, int n, int order)
{
    double zarg, zai, sinzai, infl;
    int    maxk, i;

    if (order > 8)
    {
        gmx_fatal(FARGS, "The current P3M code only supports orders up to 8");
    }

    zarg = M_PI/n;

    maxk = (n + 1)/2;

    for (i = -maxk; i < 0; i++)
    {
        zai          = zarg*i;
        sinzai       = sin(zai);
        infl         = do_p3m_influence(sinzai, order);
        bsp_mod[n+i] = infl*infl*pow(sinzai/zai, -2.0*order);
    }
    bsp_mod[0] = 1.0;
    for (i = 1; i < maxk; i++)
    {
        zai        = zarg*i;
        sinzai     = sin(zai);
        infl       = do_p3m_influence(sinzai, order);
        bsp_mod[i] = infl*infl*pow(sinzai/zai, -2.0*order);
    }
}

/* Calculate the P3M B-spline moduli */
void make_p3m_bspline_moduli(splinevec bsp_mod,
                             int nx, int ny, int nz, int order)
{
    make_p3m_bspline_moduli_dim(bsp_mod[XX], nx, order);
    make_p3m_bspline_moduli_dim(bsp_mod[YY], ny, order);
    make_p3m_bspline_moduli_dim(bsp_mod[ZZ], nz, order);
}
