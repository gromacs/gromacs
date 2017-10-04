/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "surfacearea.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

using namespace gmx;

#define UNSP_ICO_DOD      9
#define UNSP_ICO_ARC     10

#define FOURPI (4.*M_PI)
#define TORAD(A)     ((A)*0.017453293)
#define DP_TOL     0.001

static real safe_asin(real f)
{
    if ( (fabs(f) < 1.00) )
    {
        return( asin(f) );
    }
    GMX_ASSERT(fabs(f) - 1.0 > DP_TOL, "Invalid argument");
    return(M_PI_2);
}

/* routines for dot distributions on the surface of the unit sphere */
static real icosaeder_vertices(real *xus)
{
    const real rh = std::sqrt(1.-2.*cos(TORAD(72.)))/(1.-cos(TORAD(72.)));
    const real rg = cos(TORAD(72.))/(1.-cos(TORAD(72.)));
    /* icosaeder vertices */
    xus[ 0] = 0.;                  xus[ 1] = 0.;                  xus[ 2] = 1.;
    xus[ 3] = rh*cos(TORAD(72.));  xus[ 4] = rh*sin(TORAD(72.));  xus[ 5] = rg;
    xus[ 6] = rh*cos(TORAD(144.)); xus[ 7] = rh*sin(TORAD(144.)); xus[ 8] = rg;
    xus[ 9] = rh*cos(TORAD(216.)); xus[10] = rh*sin(TORAD(216.)); xus[11] = rg;
    xus[12] = rh*cos(TORAD(288.)); xus[13] = rh*sin(TORAD(288.)); xus[14] = rg;
    xus[15] = rh;                  xus[16] = 0;                   xus[17] = rg;
    xus[18] = rh*cos(TORAD(36.));  xus[19] = rh*sin(TORAD(36.));  xus[20] = -rg;
    xus[21] = rh*cos(TORAD(108.)); xus[22] = rh*sin(TORAD(108.)); xus[23] = -rg;
    xus[24] = -rh;                 xus[25] = 0;                   xus[26] = -rg;
    xus[27] = rh*cos(TORAD(252.)); xus[28] = rh*sin(TORAD(252.)); xus[29] = -rg;
    xus[30] = rh*cos(TORAD(324.)); xus[31] = rh*sin(TORAD(324.)); xus[32] = -rg;
    xus[33] = 0.;                  xus[34] = 0.;                  xus[35] = -1.;
    return rh;
}


static void divarc(real x1, real y1, real z1,
                   real x2, real y2, real z2,
                   int div1, int div2, real *xr, real *yr, real *zr)
{

    real xd, yd, zd, dd, d1, d2, s, x, y, z;
    real phi, sphi, cphi;

    xd = y1*z2-y2*z1;
    yd = z1*x2-z2*x1;
    zd = x1*y2-x2*y1;
    dd = std::sqrt(xd*xd+yd*yd+zd*zd);
    GMX_ASSERT(dd >= DP_TOL, "Rotation axis vector too short");

    d1 = x1*x1+y1*y1+z1*z1;
    d2 = x2*x2+y2*y2+z2*z2;
    GMX_ASSERT(d1 >= 0.5, "Vector 1 too short");
    GMX_ASSERT(d2 >= 0.5, "Vector 2 too short");

    phi  = safe_asin(dd/std::sqrt(d1*d2));
    phi  = phi*((real)div1)/((real)div2);
    sphi = sin(phi); cphi = cos(phi);
    s    = (x1*xd+y1*yd+z1*zd)/dd;

    x   = xd*s*(1.-cphi)/dd + x1 * cphi + (yd*z1-y1*zd)*sphi/dd;
    y   = yd*s*(1.-cphi)/dd + y1 * cphi + (zd*x1-z1*xd)*sphi/dd;
    z   = zd*s*(1.-cphi)/dd + z1 * cphi + (xd*y1-x1*yd)*sphi/dd;
    dd  = std::sqrt(x*x+y*y+z*z);
    *xr = x/dd; *yr = y/dd; *zr = z/dd;
}

/* densit...required dots per unit sphere */
static std::vector<real> ico_dot_arc(int densit)
{
    /* dot distribution on a unit sphere based on an icosaeder *
     * great circle average refining of icosahedral face       */

    int   i, j, k, tl, tl2, tn;
    real  a, d, x, y, z, x2, y2, z2, x3, y3, z3;
    real  xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
          xjk, yjk, zjk, xkj, ykj, zkj;

    /* calculate tessalation level */
    a = std::sqrt((((real) densit)-2.)/10.);
    const int  tess = (int) ceil(a);
    const int  ndot = 10*tess*tess+2;
    GMX_RELEASE_ASSERT(ndot >= densit, "Inconsistent surface dot formula");

    std::vector<real> xus(3*ndot);
    const real        rh = icosaeder_vertices(xus.data());

    if (tess > 1)
    {
        tn = 12;
        a  = rh*rh*2.*(1.-cos(TORAD(72.)));
        /* calculate tessalation of icosaeder edges */
        for (i = 0; i < 11; i++)
        {
            for (j = i+1; j < 12; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL)
                {
                    continue;
                }
                for (tl = 1; tl < tess; tl++)
                {
                    GMX_ASSERT(tn < ndot, "Inconsistent precomputed surface dot count");
                    divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                           xus[3*j], xus[1+3*j], xus[2+3*j],
                           tl, tess, &xus[3*tn], &xus[1+3*tn], &xus[2+3*tn]);
                    tn++;
                }
            }
        }
        /* calculate tessalation of icosaeder faces */
        for (i = 0; i < 10; i++)
        {
            for (j = i+1; j < 11; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL)
                {
                    continue;
                }

                for (k = j+1; k < 12; k++)
                {
                    x = xus[3*i]-xus[3*k];
                    y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(a-d) > DP_TOL)
                    {
                        continue;
                    }
                    x = xus[3*j]-xus[3*k];
                    y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(a-d) > DP_TOL)
                    {
                        continue;
                    }
                    for (tl = 1; tl < tess-1; tl++)
                    {
                        divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xji, &yji, &zji);
                        divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xki, &yki, &zki);

                        for (tl2 = 1; tl2 < tess-tl; tl2++)
                        {
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xij, &yij, &zij);
                            divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xkj, &ykj, &zkj);
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xik, &yik, &zik);
                            divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xjk, &yjk, &zjk);
                            divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl,
                                   &x, &y, &z);
                            divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2,
                                   &x2, &y2, &z2);
                            divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2,
                                   &x3, &y3, &z3);
                            GMX_ASSERT(tn < ndot,
                                       "Inconsistent precomputed surface dot count");
                            x           = x+x2+x3; y = y+y2+y3; z = z+z2+z3;
                            d           = std::sqrt(x*x+y*y+z*z);
                            xus[3*tn]   = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        } /* cycle tl2 */
                    }     /* cycle tl */
                }         /* cycle k */
            }             /* cycle j */
        }                 /* cycle i */
        GMX_ASSERT(tn == ndot, "Inconsistent precomputed surface dot count");
    }                     /* end of if (tess > 1) */

    return xus;
}                                  /* end of routine ico_dot_arc */

/* densit...required dots per unit sphere */
static std::vector<real> ico_dot_dod(int densit)
{
    /* dot distribution on a unit sphere based on an icosaeder *
     * great circle average refining of icosahedral face       */

    int   i, j, k, tl, tl2, tn, tess, j1, j2;
    real  a, d, x, y, z, x2, y2, z2, x3, y3, z3, ai_d, adod;
    real  xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
          xjk, yjk, zjk, xkj, ykj, zkj;

    /* calculate tesselation level */
    a     = std::sqrt((((real) densit)-2.)/30.);
    tess  = std::max((int) ceil(a), 1);
    const int ndot = 30*tess*tess+2;
    GMX_RELEASE_ASSERT(ndot >= densit, "Inconsistent surface dot formula");

    std::vector<real> xus(3*ndot);
    const real        rh = icosaeder_vertices(xus.data());

    tn = 12;
    /* square of the edge of an icosaeder */
    a = rh*rh*2.*(1.-cos(TORAD(72.)));
    /* dodecaeder vertices */
    for (i = 0; i < 10; i++)
    {
        for (j = i+1; j < 11; j++)
        {
            x = xus[3*i]-xus[3*j];
            y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
            d = x*x+y*y+z*z;
            if (fabs(a-d) > DP_TOL)
            {
                continue;
            }
            for (k = j+1; k < 12; k++)
            {
                x = xus[3*i]-xus[3*k];
                y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL)
                {
                    continue;
                }
                x = xus[3*j]-xus[3*k];
                y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL)
                {
                    continue;
                }
                x         = xus[  3*i]+xus[  3*j]+xus[  3*k];
                y         = xus[1+3*i]+xus[1+3*j]+xus[1+3*k];
                z         = xus[2+3*i]+xus[2+3*j]+xus[2+3*k];
                d         = std::sqrt(x*x+y*y+z*z);
                xus[3*tn] = x/d; xus[1+3*tn] = y/d; xus[2+3*tn] = z/d;
                tn++;
            }
        }
    }

    if (tess > 1)
    {
        tn = 32;
        /* square of the edge of an dodecaeder */
        adod = 4.*(cos(TORAD(108.))-cos(TORAD(120.)))/(1.-cos(TORAD(120.)));
        /* square of the distance of two adjacent vertices of ico- and dodecaeder */
        ai_d = 2.*(1.-std::sqrt(1.-a/3.));

        /* calculate tessalation of mixed edges */
        for (i = 0; i < 31; i++)
        {
            j1 = 12; j2 = 32; a = ai_d;
            if (i >= 12)
            {
                j1 = i+1; a = adod;
            }
            for (j = j1; j < j2; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL)
                {
                    continue;
                }
                for (tl = 1; tl < tess; tl++)
                {
                    GMX_ASSERT(tn < ndot,
                               "Inconsistent precomputed surface dot count");
                    divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                           xus[3*j], xus[1+3*j], xus[2+3*j],
                           tl, tess, &xus[3*tn], &xus[1+3*tn], &xus[2+3*tn]);
                    tn++;
                }
            }
        }
        /* calculate tessalation of pentakisdodecahedron faces */
        for (i = 0; i < 12; i++)
        {
            for (j = 12; j < 31; j++)
            {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(ai_d-d) > DP_TOL)
                {
                    continue;
                }

                for (k = j+1; k < 32; k++)
                {
                    x = xus[3*i]-xus[3*k];
                    y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(ai_d-d) > DP_TOL)
                    {
                        continue;
                    }
                    x = xus[3*j]-xus[3*k];
                    y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(adod-d) > DP_TOL)
                    {
                        continue;
                    }
                    for (tl = 1; tl < tess-1; tl++)
                    {
                        divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xji, &yji, &zji);
                        divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                               xus[3*i], xus[1+3*i], xus[2+3*i],
                               tl, tess, &xki, &yki, &zki);

                        for (tl2 = 1; tl2 < tess-tl; tl2++)
                        {
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xij, &yij, &zij);
                            divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                                   xus[3*j], xus[1+3*j], xus[2+3*j],
                                   tl2, tess, &xkj, &ykj, &zkj);
                            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xik, &yik, &zik);
                            divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                                   xus[3*k], xus[1+3*k], xus[2+3*k],
                                   tess-tl-tl2, tess, &xjk, &yjk, &zjk);
                            divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl,
                                   &x, &y, &z);
                            divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2,
                                   &x2, &y2, &z2);
                            divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2,
                                   &x3, &y3, &z3);
                            GMX_ASSERT(tn < ndot,
                                       "Inconsistent precomputed surface dot count");
                            x           = x+x2+x3; y = y+y2+y3; z = z+z2+z3;
                            d           = std::sqrt(x*x+y*y+z*z);
                            xus[3*tn]   = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        } /* cycle tl2 */
                    }     /* cycle tl */
                }         /* cycle k */
            }             /* cycle j */
        }                 /* cycle i */
        GMX_ASSERT(tn == ndot, "Inconsistent precomputed surface dot count");
    }                     /* end of if (tess > 1) */

    return xus;
}       /* end of routine ico_dot_dod */

static int unsp_type(int densit)
{
    int i1, i2;
    i1 = 1;
    while (10*i1*i1+2 < densit)
    {
        i1++;
    }
    i2 = 1;
    while (30*i2*i2+2 < densit)
    {
        i2++;
    }
    if (10*i1*i1-2 < 30*i2*i2-2)
    {
        return UNSP_ICO_ARC;
    }
    else
    {
        return UNSP_ICO_DOD;
    }
}

static std::vector<real> make_unsp(int densit, int cubus)
{
    int              *ico_wk, *ico_pt;
    int               ico_cube, ico_cube_cb, i, j, k, l, ijk, tn, tl, tl2;
    int              *work;
    real              x, y, z;

    int               mode = unsp_type(densit);
    std::vector<real> xus;
    if (mode == UNSP_ICO_ARC)
    {
        xus = ico_dot_arc(densit);
    }
    else if (mode == UNSP_ICO_DOD)
    {
        xus = ico_dot_dod(densit);
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Invalid unit sphere mode");
    }

    const int ndot = static_cast<int>(xus.size())/3;

    /* determine distribution of points in elementary cubes */
    if (cubus)
    {
        ico_cube = cubus;
    }
    else
    {
        i          = 1;
        while (i*i*i*2 < ndot)
        {
            i++;
        }
        ico_cube = std::max(i-1, 0);
    }
    ico_cube_cb = ico_cube*ico_cube*ico_cube;
    const real del_cube = 2./((real)ico_cube);
    snew(work, ndot);
    for (l = 0; l < ndot; l++)
    {
        i = std::max((int) floor((1.+xus[3*l])/del_cube), 0);
        if (i >= ico_cube)
        {
            i = ico_cube-1;
        }
        j = std::max((int) floor((1.+xus[1+3*l])/del_cube), 0);
        if (j >= ico_cube)
        {
            j = ico_cube-1;
        }
        k = std::max((int) floor((1.+xus[2+3*l])/del_cube), 0);
        if (k >= ico_cube)
        {
            k = ico_cube-1;
        }
        ijk     = i+j*ico_cube+k*ico_cube*ico_cube;
        work[l] = ijk;
    }

    snew(ico_wk, 2*ico_cube_cb+1);

    ico_pt = ico_wk+ico_cube_cb;
    for (l = 0; l < ndot; l++)
    {
        ico_wk[work[l]]++; /* dots per elementary cube */
    }

    /* reordering of the coordinate array in accordance with box number */
    tn = 0;
    for (i = 0; i < ico_cube; i++)
    {
        for (j = 0; j < ico_cube; j++)
        {
            for (k = 0; k < ico_cube; k++)
            {
                tl            = 0;
                tl2           = tn;
                ijk           = i+ico_cube*j+ico_cube*ico_cube*k;
                *(ico_pt+ijk) = tn;
                for (l = tl2; l < ndot; l++)
                {
                    if (ijk == work[l])
                    {
                        x          = xus[3*l]; y = xus[1+3*l]; z = xus[2+3*l];
                        xus[3*l]   = xus[3*tn];
                        xus[1+3*l] = xus[1+3*tn]; xus[2+3*l] = xus[2+3*tn];
                        xus[3*tn]  = x; xus[1+3*tn] = y; xus[2+3*tn] = z;
                        ijk        = work[l]; work[l] = work[tn]; work[tn] = ijk;
                        tn++; tl++;
                    }
                }
                *(ico_wk+ijk) = tl;
            } /* cycle k */
        }     /* cycle j */
    }         /* cycle i */
    sfree(ico_wk);
    sfree(work);

    return xus;
}

static void
nsc_dclm_pbc(const rvec *coords, const ArrayRef<const real> &radius, int nat,
             const real *xus, int n_dot, int mode,
             real *value_of_area, real **at_area,
             real *value_of_vol,
             real **lidots, int *nu_dots,
             int index[], AnalysisNeighborhood *nb,
             const t_pbc *pbc)
{
    const real dotarea = FOURPI/(real) n_dot;

    if (debug)
    {
        fprintf(debug, "nsc_dclm: n_dot=%5d %9.3f\n", n_dot, dotarea);
    }

    /* start with neighbour list */
    /* calculate neighbour list with the box algorithm */
    if (nat == 0)
    {
        return;
    }
    real        area = 0.0, vol = 0.0;
    real       *dots = nullptr, *atom_area = nullptr;
    int         lfnr = 0, maxdots = 0;
    if (mode & FLAG_VOLUME)
    {
        vol = 0.;
    }
    if (mode & FLAG_DOTS)
    {
        maxdots = (3*n_dot*nat)/10;
        snew(dots, maxdots);
        lfnr = 0;
    }
    if (mode & FLAG_ATOM_AREA)
    {
        snew(atom_area, nat);
    }

    // Compute the center of the molecule for volume calculation.
    // In principle, the center should not influence the results, but that is
    // only true at the limit of infinite dot density, so this makes the
    // results translation-invariant.
    // With PBC, if the molecule is broken across the boundary, the computation
    // is broken in other ways as well, so it does not need to be considered
    // here.
    real xs = 0.0, ys = 0.0, zs = 0.0;
    for (int i = 0; i < nat; ++i)
    {
        const int iat = index[i];
        xs += coords[iat][XX];
        ys += coords[iat][YY];
        zs += coords[iat][ZZ];
    }
    xs /= nat;
    ys /= nat;
    zs /= nat;

    AnalysisNeighborhoodPositions pos(coords, radius.size());
    pos.indexed(constArrayRefFromArray(index, nat));
    AnalysisNeighborhoodSearch    nbsearch(nb->initSearch(pbc, pos));

    std::vector<int>              wkdot(n_dot);

    for (int i = 0; i < nat; ++i)
    {
        const int                      iat  = index[i];
        const real                     ai   = radius[iat];
        const real                     aisq = ai*ai;
        AnalysisNeighborhoodPairSearch pairSearch(
                nbsearch.startPairSearch(coords[iat]));
        AnalysisNeighborhoodPair       pair;
        std::fill(wkdot.begin(), wkdot.end(), 1);
        int currDotCount = n_dot;
        while (currDotCount > 0 && pairSearch.findNextPair(&pair))
        {
            const int  jat = index[pair.refIndex()];
            const real aj  = radius[jat];
            const real d2  = pair.distance2();
            if (iat == jat || d2 > gmx::square(ai+aj))
            {
                continue;
            }
            const rvec &dx     = pair.dx();
            const real  refdot = (d2 + aisq - aj*aj)/(2*ai);
            // TODO: Consider whether micro-optimizations from the old
            // implementation would be useful, compared to the complexity that
            // they bring: instead of this direct loop, the neighbors were
            // stored into a temporary array, the loop order was
            // reversed (first over dots, then over neighbors), and for each
            // dot, it was first checked whether the same neighbor that
            // resulted in marking the previous dot covered would also cover
            // this dot. This presumably plays together with sorting of the
            // surface dots (done in make_unsp) to avoid some of the looping.
            // Alternatively, we could keep a skip list here to avoid
            // repeatedly looping over dots that have already marked as
            // covered.
            for (int j = 0; j < n_dot; ++j)
            {
                if (wkdot[j] && iprod(&xus[3*j], dx) > refdot)
                {
                    --currDotCount;
                    wkdot[j] = 0;
                }
            }
        }

        const real a = aisq * dotarea * currDotCount;
        area = area + a;
        if (mode & FLAG_ATOM_AREA)
        {
            atom_area[i] = a;
        }
        const real xi = coords[iat][XX];
        const real yi = coords[iat][YY];
        const real zi = coords[iat][ZZ];
        if (mode & FLAG_DOTS)
        {
            for (int l = 0; l < n_dot; l++)
            {
                if (wkdot[l])
                {
                    lfnr++;
                    if (maxdots <= 3*lfnr+1)
                    {
                        maxdots = maxdots+n_dot*3;
                        srenew(dots, maxdots);
                    }
                    dots[3*lfnr-3] = ai*xus[3*l]+xi;
                    dots[3*lfnr-2] = ai*xus[1+3*l]+yi;
                    dots[3*lfnr-1] = ai*xus[2+3*l]+zi;
                }
            }
        }
        if (mode & FLAG_VOLUME)
        {
            real dx = 0.0, dy = 0.0, dz = 0.0;
            for (int l = 0; l < n_dot; l++)
            {
                if (wkdot[l])
                {
                    dx = dx+xus[3*l];
                    dy = dy+xus[1+3*l];
                    dz = dz+xus[2+3*l];
                }
            }
            vol = vol+aisq*(dx*(xi-xs)+dy*(yi-ys)+dz*(zi-zs) + ai*currDotCount);
        }
    }

    if (mode & FLAG_VOLUME)
    {
        *value_of_vol = vol*FOURPI/(3.*n_dot);
    }
    if (mode & FLAG_DOTS)
    {
        *nu_dots = lfnr;
        *lidots  = dots;
    }
    if (mode & FLAG_ATOM_AREA)
    {
        *at_area = atom_area;
    }
    *value_of_area = area;

    if (debug)
    {
        fprintf(debug, "area=%8.3f\n", area);
    }
}

namespace gmx
{

class SurfaceAreaCalculator::Impl
{
    public:
        Impl() : flags_(0)
        {
        }

        std::vector<real>             unitSphereDots_;
        ArrayRef<const real>          radius_;
        int                           flags_;
        mutable AnalysisNeighborhood  nb_;
};

SurfaceAreaCalculator::SurfaceAreaCalculator()
    : impl_(new Impl())
{
}

SurfaceAreaCalculator::~SurfaceAreaCalculator()
{
}

void SurfaceAreaCalculator::setDotCount(int dotCount)
{
    impl_->unitSphereDots_ = make_unsp(dotCount, 4);
}

void SurfaceAreaCalculator::setRadii(const ArrayRef<const real> &radius)
{
    impl_->radius_ = radius;
    if (!radius.empty())
    {
        const real maxRadius = *std::max_element(radius.begin(), radius.end());
        impl_->nb_.setCutoff(2*maxRadius);
    }
}

void SurfaceAreaCalculator::setCalculateVolume(bool bVolume)
{
    if (bVolume)
    {
        impl_->flags_ |= FLAG_VOLUME;
    }
    else
    {
        impl_->flags_ &= ~FLAG_VOLUME;
    }
}

void SurfaceAreaCalculator::setCalculateAtomArea(bool bAtomArea)
{
    if (bAtomArea)
    {
        impl_->flags_ |= FLAG_ATOM_AREA;
    }
    else
    {
        impl_->flags_ &= ~FLAG_ATOM_AREA;
    }
}

void SurfaceAreaCalculator::setCalculateSurfaceDots(bool bDots)
{
    if (bDots)
    {
        impl_->flags_ |= FLAG_DOTS;
    }
    else
    {
        impl_->flags_ &= ~FLAG_DOTS;
    }
}

void SurfaceAreaCalculator::calculate(
        const rvec *x, const t_pbc *pbc,
        int nat, int index[], int flags, real *area, real *volume,
        real **at_area, real **lidots, int *n_dots) const
{
    flags |= impl_->flags_;
    *area  = 0;
    if (volume == nullptr)
    {
        flags &= ~FLAG_VOLUME;
    }
    else
    {
        *volume = 0;
    }
    if (at_area == nullptr)
    {
        flags &= ~FLAG_ATOM_AREA;
    }
    else
    {
        *at_area = nullptr;
    }
    if (lidots == nullptr)
    {
        flags &= ~FLAG_DOTS;
    }
    else
    {
        *lidots = nullptr;
    }
    if (n_dots == nullptr)
    {
        flags &= ~FLAG_DOTS;
    }
    else
    {
        *n_dots = 0;
    }
    nsc_dclm_pbc(x, impl_->radius_, nat,
                 &impl_->unitSphereDots_[0], impl_->unitSphereDots_.size()/3,
                 flags, area, at_area, volume, lidots, n_dots, index,
                 &impl_->nb_, pbc);
}

} // namespace gmx
