/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
/* Modified DvdS */
#include "pbc.h"
#include "macros.h"
#include "vec.h"
#include "smalloc.h"
#include "nsc.h"

#define TEST_NSC 0

#define TEST_ARC 0
#define TEST_DOD 0
#define TEST_CUBE 0

#define UNSP_ICO_DOD      9
#define UNSP_ICO_ARC     10

real   *xpunsp = NULL;
real    del_cube;
int    *ico_wk = NULL, *ico_pt = NULL;
int     n_dot, ico_cube, last_n_dot = 0, last_densit = 0, last_unsp = 0;
int     last_cubus = 0;

#define FOURPI (4.*M_PI)
#define TORAD(A)     ((A)*0.017453293)
#define DP_TOL     0.001

#define UPDATE_FL  __file__ = __FILE__, __line__ = __LINE__
const char * __file__;   /* declared versions of macros */
int          __line__;   /* __FILE__  and __LINE__ */

#ifdef ERROR
#undef ERROR
#endif
#define ERROR UPDATE_FL, error
void error(const char *fmt, ...)
{
    va_list args;
    fprintf(stderr,
            "\n---> ERROR when executing line %i in file %s !\n",
            __line__, __file__);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n---> Execution stopped !\n\n");
}

#define WARNING UPDATE_FL, warning2
void warning2(const char *fmt, ...)
{
    va_list args;
    fprintf(stderr,
            "\n---> WARNING : line %i in file %s\n",
            __line__, __file__);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, " ...!\n\n");
    fflush(stderr);
    fflush(stdout);
}

#define ASIN safe_asin
real safe_asin(real f)
{
    if ( (fabs(f) < 1.00) )
    {
        return( asin(f) );
    }
    if ( (fabs(f) - 1.00)  <= DP_TOL)
    {
        ERROR("ASIN : invalid argument %f", f);
    }
    return(M_PI_2);
}




/* routines for dot distributions on the surface of the unit sphere */
real rg, rh;

void icosaeder_vertices(real *xus)
{
    rh = sqrt(1.-2.*cos(TORAD(72.)))/(1.-cos(TORAD(72.)));
    rg = cos(TORAD(72.))/(1.-cos(TORAD(72.)));
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
}


void divarc(real x1, real y1, real z1,
            real x2, real y2, real z2,
            int div1, int div2, real *xr, real *yr, real *zr)
{

    real xd, yd, zd, dd, d1, d2, s, x, y, z;
    real phi, sphi, cphi;

    xd = y1*z2-y2*z1;
    yd = z1*x2-z2*x1;
    zd = x1*y2-x2*y1;
    dd = sqrt(xd*xd+yd*yd+zd*zd);
    if (dd < DP_TOL)
    {
        ERROR("divarc: rotation axis of length %f", dd);
    }

    d1 = x1*x1+y1*y1+z1*z1;
    if (d1 < 0.5)
    {
        ERROR("divarc: vector 1 of sq.length %f", d1);
    }
    d2 = x2*x2+y2*y2+z2*z2;
    if (d2 < 0.5)
    {
        ERROR("divarc: vector 2 of sq.length %f", d2);
    }

    phi  = ASIN(dd/sqrt(d1*d2));
    phi  = phi*((real)div1)/((real)div2);
    sphi = sin(phi); cphi = cos(phi);
    s    = (x1*xd+y1*yd+z1*zd)/dd;

    x   = xd*s*(1.-cphi)/dd + x1 * cphi + (yd*z1-y1*zd)*sphi/dd;
    y   = yd*s*(1.-cphi)/dd + y1 * cphi + (zd*x1-z1*xd)*sphi/dd;
    z   = zd*s*(1.-cphi)/dd + z1 * cphi + (xd*y1-x1*yd)*sphi/dd;
    dd  = sqrt(x*x+y*y+z*z);
    *xr = x/dd; *yr = y/dd; *zr = z/dd;
}

int ico_dot_arc(int densit)   /* densit...required dots per unit sphere */
{
    /* dot distribution on a unit sphere based on an icosaeder *
     * great circle average refining of icosahedral face       */

    int   i, j, k, tl, tl2, tn, tess;
    real  a, d, x, y, z, x2, y2, z2, x3, y3, z3;
    real  xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
          xjk, yjk, zjk, xkj, ykj, zkj;
    real *xus = NULL;

    /* calculate tessalation level */
    a     = sqrt((((real) densit)-2.)/10.);
    tess  = (int) ceil(a);
    n_dot = 10*tess*tess+2;
    if (n_dot < densit)
    {
        ERROR("ico_dot_arc: error in formula for tessalation level (%d->%d, %d)",
              tess, n_dot, densit);
    }

    snew(xus, 3*n_dot);
    xpunsp = xus;
    icosaeder_vertices(xus);

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
                    if (tn >= n_dot)
                    {
                        ERROR("ico_dot: tn exceeds dimension of xus");
                    }
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
                            if (tn >= n_dot)
                            {
                                ERROR("ico_dot: tn exceeds dimension of xus");
                            }
                            divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl,
                                   &x, &y, &z);
                            divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2,
                                   &x2, &y2, &z2);
                            divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2,
                                   &x3, &y3, &z3);
                            x           = x+x2+x3; y = y+y2+y3; z = z+z2+z3;
                            d           = sqrt(x*x+y*y+z*z);
                            xus[3*tn]   = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        } /* cycle tl2 */
                    }     /* cycle tl */
                }         /* cycle k */
            }             /* cycle j */
        }                 /* cycle i */
        if (n_dot != tn)
        {
            ERROR("ico_dot: n_dot(%d) and tn(%d) differ", n_dot, tn);
        }
    }   /* end of if (tess > 1) */

    return n_dot;
}                           /* end of routine ico_dot_arc */

int ico_dot_dod(int densit) /* densit...required dots per unit sphere */
{
    /* dot distribution on a unit sphere based on an icosaeder *
     * great circle average refining of icosahedral face       */

    int   i, j, k, tl, tl2, tn, tess, j1, j2;
    real  a, d, x, y, z, x2, y2, z2, x3, y3, z3, ai_d, adod;
    real  xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
          xjk, yjk, zjk, xkj, ykj, zkj;
    real *xus = NULL;
    /* calculate tesselation level */
    a     = sqrt((((real) densit)-2.)/30.);
    tess  = max((int) ceil(a), 1);
    n_dot = 30*tess*tess+2;
    if (n_dot < densit)
    {
        ERROR("ico_dot_dod: error in formula for tessalation level (%d->%d, %d)",
              tess, n_dot, densit);
    }

    snew(xus, 3*n_dot);
    xpunsp = xus;
    icosaeder_vertices(xus);

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
                d         = sqrt(x*x+y*y+z*z);
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
        ai_d = 2.*(1.-sqrt(1.-a/3.));

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
                    if (tn >= n_dot)
                    {
                        ERROR("ico_dot: tn exceeds dimension of xus");
                    }
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
                            if (tn >= n_dot)
                            {
                                ERROR("ico_dot: tn exceeds dimension of xus");
                            }
                            divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl,
                                   &x, &y, &z);
                            divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2,
                                   &x2, &y2, &z2);
                            divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2,
                                   &x3, &y3, &z3);
                            x           = x+x2+x3; y = y+y2+y3; z = z+z2+z3;
                            d           = sqrt(x*x+y*y+z*z);
                            xus[3*tn]   = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        } /* cycle tl2 */
                    }     /* cycle tl */
                }         /* cycle k */
            }             /* cycle j */
        }                 /* cycle i */
        if (n_dot != tn)
        {
            ERROR("ico_dot: n_dot(%d) and tn(%d) differ", n_dot, tn);
        }
    }   /* end of if (tess > 1) */

    return n_dot;
}       /* end of routine ico_dot_dod */

int unsp_type(int densit)
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

int make_unsp(int densit, int mode, int * num_dot, int cubus)
{
    int   ndot, ico_cube_cb, i, j, k, l, ijk, tn, tl, tl2;
    real *xus;
    int  *work;
    real  x, y, z;

    if (xpunsp)
    {
        free(xpunsp);
    }
    if (ico_wk)
    {
        free(ico_wk);
    }

    k = 1; if (mode < 0)
    {
        k = 0; mode = -mode;
    }
    if (mode == UNSP_ICO_ARC)
    {
        ndot = ico_dot_arc(densit);
    }
    else if (mode == UNSP_ICO_DOD)
    {
        ndot = ico_dot_dod(densit);
    }
    else
    {
        WARNING("make_unsp: mode %c%d not allowed", (k) ? '+' : '-', mode);
        return 1;
    }

    last_n_dot = ndot; last_densit = densit; last_unsp = mode;
    *num_dot   = ndot; if (k)
    {
        return 0;
    }

    /* in the following the dots of the unit sphere may be resorted */
    last_unsp = -last_unsp;

    /* determine distribution of points in elementary cubes */
    if (cubus)
    {
        ico_cube = cubus;
    }
    else
    {
        last_cubus = 0;
        i          = 1;
        while (i*i*i*2 < ndot)
        {
            i++;
        }
        ico_cube = max(i-1, 0);
    }
    ico_cube_cb = ico_cube*ico_cube*ico_cube;
    del_cube    = 2./((real)ico_cube);
    snew(work, ndot);
    xus = xpunsp;
    for (l = 0; l < ndot; l++)
    {
        i = max((int) floor((1.+xus[3*l])/del_cube), 0);
        if (i >= ico_cube)
        {
            i = ico_cube-1;
        }
        j = max((int) floor((1.+xus[1+3*l])/del_cube), 0);
        if (j >= ico_cube)
        {
            j = ico_cube-1;
        }
        k = max((int) floor((1.+xus[2+3*l])/del_cube), 0);
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
    free(work);

    return 0;
}


typedef struct _stwknb {
    real x;
    real y;
    real z;
    real dot;
} Neighb;

int nsc_dclm_pbc(rvec *coords, real *radius, int nat,
                 int  densit, int mode,
                 real *value_of_area, real **at_area,
                 real *value_of_vol,
                 real **lidots, int *nu_dots,
                 atom_id index[], int ePBC, matrix box)
{
    int      iat, i, ii, iii, ix, iy, iz, ixe, ixs, iye, iys, ize, izs, i_ac;
    int      jat, j, jj, jjj, jx, jy, jz;
    int      distribution;
    int      l;
    int      maxnei, nnei, last, maxdots = 0;
    int     *wkdot = NULL, *wkbox = NULL, *wkat1 = NULL, *wkatm = NULL;
    Neighb  *wknb, *ctnb;
    int      iii1, iii2, iiat, lfnr = 0, i_at, j_at;
    real     dx, dy, dz, dd, ai, aisq, ajsq, aj, as, a;
    real     xi, yi, zi, xs = 0., ys = 0., zs = 0.;
    real     dotarea, area, vol = 0.;
    real    *xus, *dots = NULL, *atom_area = NULL;
    int      nxbox, nybox, nzbox, nxy, nxyz;
    real     xmin = 0, ymin = 0, zmin = 0, xmax, ymax, zmax, ra2max, d, *pco;
    /* Added DvdS 2006-07-19 */
    t_pbc    pbc;
    rvec     ddx, *x = NULL;
    int      iat_xx, jat_xx;

    distribution = unsp_type(densit);
    if (distribution != -last_unsp || last_cubus != 4 ||
        (densit != last_densit && densit != last_n_dot))
    {
        if (make_unsp(densit, (-distribution), &n_dot, 4))
        {
            return 1;
        }
    }
    xus = xpunsp;

    dotarea = FOURPI/(real) n_dot;
    area    = 0.;

    if (debug)
    {
        fprintf(debug, "nsc_dclm: n_dot=%5d %9.3f\n", n_dot, dotarea);
    }

    /* start with neighbour list */
    /* calculate neighbour list with the box algorithm */
    if (nat == 0)
    {
        WARNING("nsc_dclm: no surface atoms selected");
        return 1;
    }
    if (mode & FLAG_VOLUME)
    {
        vol = 0.;
    }
    if (mode & FLAG_DOTS)
    {
        maxdots = (3*n_dot*nat)/10;
        /* should be set to NULL on first user call */
        if (dots == NULL)
        {
            snew(dots, maxdots);
        }
        else
        {
            srenew(dots, maxdots);
        }

        lfnr = 0;
    }
    if (mode & FLAG_ATOM_AREA)
    {
        /* should be set to NULL on first user call */
        if (atom_area == NULL)
        {
            snew(atom_area, nat);
        }
        else
        {
            srenew(atom_area, nat);
        }
    }

    /* Compute minimum size for grid cells */
    ra2max = radius[index[0]];
    for (iat_xx = 1; (iat_xx < nat); iat_xx++)
    {
        iat    = index[iat_xx];
        ra2max = max(ra2max, radius[iat]);
    }
    ra2max = 2*ra2max;

    /* Added DvdS 2006-07-19 */
    /* Updated 2008-10-09 */
    if (box)
    {
        set_pbc(&pbc, ePBC, box);
        snew(x, nat);
        for (i = 0; (i < nat); i++)
        {
            iat  = index[0];
            copy_rvec(coords[iat], x[i]);
        }
        put_atoms_in_triclinic_unitcell(ecenterTRIC, box, nat, x);
        nxbox = max(1, floor(norm(box[XX])/ra2max));
        nybox = max(1, floor(norm(box[YY])/ra2max));
        nzbox = max(1, floor(norm(box[ZZ])/ra2max));
        if (debug)
        {
            fprintf(debug, "nbox = %d, %d, %d\n", nxbox, nybox, nzbox);
        }
    }
    else
    {
        /* dimensions of atomic set, cell edge is 2*ra_max */
        iat    = index[0];
        xmin   = coords[iat][XX]; xmax = xmin; xs = xmin;
        ymin   = coords[iat][YY]; ymax = ymin; ys = ymin;
        zmin   = coords[iat][ZZ]; zmax = zmin; zs = zmin;

        for (iat_xx = 1; (iat_xx < nat); iat_xx++)
        {
            iat  = index[iat_xx];
            pco  = coords[iat];
            xmin = min(xmin, *pco);     xmax = max(xmax, *pco);
            ymin = min(ymin, *(pco+1)); ymax = max(ymax, *(pco+1));
            zmin = min(zmin, *(pco+2)); zmax = max(zmax, *(pco+2));
            xs   = xs+ *pco; ys = ys+ *(pco+1); zs = zs+ *(pco+2);
        }
        xs = xs/ (real) nat;
        ys = ys/ (real) nat;
        zs = zs/ (real) nat;
        if (debug)
        {
            fprintf(debug, "nsc_dclm: n_dot=%5d ra2max=%9.3f %9.3f\n", n_dot, ra2max, dotarea);
        }

        d    = xmax-xmin; nxbox = (int) max(ceil(d/ra2max), 1.);
        d    = (((real)nxbox)*ra2max-d)/2.;
        xmin = xmin-d; xmax = xmax+d;
        d    = ymax-ymin; nybox = (int) max(ceil(d/ra2max), 1.);
        d    = (((real)nybox)*ra2max-d)/2.;
        ymin = ymin-d; ymax = ymax+d;
        d    = zmax-zmin; nzbox = (int) max(ceil(d/ra2max), 1.);
        d    = (((real)nzbox)*ra2max-d)/2.;
        zmin = zmin-d; zmax = zmax+d;
    }
    /* Help variables */
    nxy  = nxbox*nybox;
    nxyz = nxy*nzbox;

    /* box number of atoms */
    snew(wkatm, nat);
    snew(wkat1, nat);
    snew(wkdot, n_dot);
    snew(wkbox, nxyz+1);

    if (box)
    {
        matrix box_1;
        rvec   x_1;
        int    ix, iy, iz, m;
        m_inv(box, box_1);
        for (i = 0; (i < nat); i++)
        {
            mvmul(box_1, x[i], x_1);
            ix = ((int)floor(x_1[XX]*nxbox) + 2*nxbox) % nxbox;
            iy = ((int)floor(x_1[YY]*nybox) + 2*nybox) % nybox;
            iz = ((int)floor(x_1[ZZ]*nzbox) + 2*nzbox) % nzbox;
            j  =  ix + iy*nxbox + iz*nxbox*nybox;
            if (debug)
            {
                fprintf(debug, "Atom %d cell index %d. x = (%8.3f,%8.3f,%8.3f) fc = (%8.3f,%8.3f,%8.3f)\n",
                        i, j, x[i][XX], x[i][YY], x[i][ZZ], x_1[XX], x_1[YY], x_1[ZZ]);
            }
            range_check(j, 0, nxyz);
            wkat1[i] = j;
            wkbox[j]++;
        }
    }
    else
    {
        /* Put the atoms in their boxes */
        for (iat_xx = 0; (iat_xx < nat); iat_xx++)
        {
            iat           = index[iat_xx];
            pco           = coords[iat];
            i             = (int) max(floor((pco[XX]-xmin)/ra2max), 0);
            i             = min(i, nxbox-1);
            j             = (int) max(floor((pco[YY]-ymin)/ra2max), 0);
            j             = min(j, nybox-1);
            l             = (int) max(floor((pco[ZZ]-zmin)/ra2max), 0);
            l             = min(l, nzbox-1);
            i             = i+j*nxbox+l*nxy;
            wkat1[iat_xx] = i;
            wkbox[i]++;
        }
    }

    /* sorting of atoms in accordance with box numbers */
    j = wkbox[0];
    for (i = 1; i < nxyz; i++)
    {
        j = max(wkbox[i], j);
    }
    for (i = 1; i <= nxyz; i++)
    {
        wkbox[i] += wkbox[i-1];
    }

    /* maxnei = (int) floor(ra2max*ra2max*ra2max*0.5); */
    maxnei = min(nat, 27*j);
    snew(wknb, maxnei);
    for (iat_xx = 0; iat_xx < nat; iat_xx++)
    {
        iat = index[iat_xx];
        range_check(wkat1[iat_xx], 0, nxyz);
        wkatm[--wkbox[wkat1[iat_xx]]] = iat_xx;
        if (debug)
        {
            fprintf(debug, "atom %5d on place %5d\n", iat, wkbox[wkat1[iat_xx]]);
        }
    }

    if (debug)
    {
        fprintf(debug, "nsc_dclm: n_dot=%5d ra2max=%9.3f %9.3f\n",
                n_dot, ra2max, dotarea);
        fprintf(debug, "neighbour list calculated/box(xyz):%d %d %d\n",
                nxbox, nybox, nzbox);

        for (i = 0; i < nxyz; i++)
        {
            fprintf(debug, "box %6d : atoms %4d-%4d    %5d\n",
                    i, wkbox[i], wkbox[i+1]-1, wkbox[i+1]-wkbox[i]);
        }
        for (i = 0; i < nat; i++)
        {
            fprintf(debug, "list place %5d by atom %7d\n", i, index[wkatm[i]]);
        }
    }

    /* calculate surface for all atoms, step cube-wise */
    for (iz = 0; iz < nzbox; iz++)
    {
        iii = iz*nxy;
        if (box)
        {
            izs = iz-1;
            ize = min(iz+2, izs+nzbox);
        }
        else
        {
            izs = max(iz-1, 0);
            ize = min(iz+2, nzbox);
        }
        for (iy = 0; iy < nybox; iy++)
        {
            ii = iy*nxbox+iii;
            if (box)
            {
                iys = iy-1;
                iye = min(iy+2, iys+nybox);
            }
            else
            {
                iys = max(iy-1, 0);
                iye = min(iy+2, nybox);
            }
            for (ix = 0; ix < nxbox; ix++)
            {
                i    = ii+ix;
                iii1 = wkbox[i];
                iii2 = wkbox[i+1];
                if (iii1 >= iii2)
                {
                    continue;
                }
                if (box)
                {
                    ixs = ix-1;
                    ixe = min(ix+2, ixs+nxbox);
                }
                else
                {
                    ixs = max(ix-1, 0);
                    ixe = min(ix+2, nxbox);
                }
                iiat = 0;
                /* make intermediate atom list */
                for (jz = izs; jz < ize; jz++)
                {
                    jjj = ((jz+nzbox) % nzbox)*nxy;
                    for (jy = iys; jy < iye; jy++)
                    {
                        jj = ((jy+nybox) % nybox)*nxbox+jjj;
                        for (jx = ixs; jx < ixe; jx++)
                        {
                            j = jj+((jx+nxbox) % nxbox);
                            for (jat = wkbox[j]; jat < wkbox[j+1]; jat++)
                            {
                                range_check(wkatm[jat], 0, nat);
                                range_check(iiat, 0, nat);
                                wkat1[iiat] = wkatm[jat];
                                iiat++;
                            } /* end of cycle "jat" */
                        }     /* end of cycle "jx" */
                    }         /* end of cycle "jy" */
                }             /* end of cycle "jz" */
                for (iat = iii1; iat < iii2; iat++)
                {
                    i_at = index[wkatm[iat]];
                    ai   = radius[i_at];
                    aisq = ai*ai;
                    pco  = coords[i_at];
                    xi   = pco[XX]; yi = pco[YY]; zi = pco[ZZ];
                    for (i = 0; i < n_dot; i++)
                    {
                        wkdot[i] = 0;
                    }

                    ctnb = wknb; nnei = 0;
                    for (j = 0; j < iiat; j++)
                    {
                        j_at = index[wkat1[j]];
                        if (j_at == i_at)
                        {
                            continue;
                        }

                        aj   = radius[j_at];
                        ajsq = aj*aj;
                        pco  = coords[j_at];

                        /* Added DvdS 2006-07-19 */
                        if (box)
                        {
                            /*rvec xxi;

                               xxi[XX] = xi;
                               xxi[YY] = yi;
                               xxi[ZZ] = zi;
                               pbc_dx(&pbc,pco,xxi,ddx);*/
                            pbc_dx(&pbc, coords[j_at], coords[i_at], ddx);
                            dx = ddx[XX];
                            dy = ddx[YY];
                            dz = ddx[ZZ];
                        }
                        else
                        {
                            dx = pco[XX]-xi;
                            dy = pco[YY]-yi;
                            dz = pco[ZZ]-zi;
                        }
                        dd = dx*dx+dy*dy+dz*dz;
                        as = ai+aj;
                        if (dd > as*as)
                        {
                            continue;
                        }
                        nnei++;
                        ctnb->x   = dx;
                        ctnb->y   = dy;
                        ctnb->z   = dz;
                        ctnb->dot = (dd+aisq-ajsq)/(2.*ai); /* reference dot product */
                        ctnb++;
                    }

                    /* check points on accessibility */
                    if (nnei)
                    {
                        last = 0; i_ac = 0;
                        for (l = 0; l < n_dot; l++)
                        {
                            if (xus[3*l]*(wknb+last)->x+
                                xus[1+3*l]*(wknb+last)->y+
                                xus[2+3*l]*(wknb+last)->z <= (wknb+last)->dot)
                            {
                                for (j = 0; j < nnei; j++)
                                {
                                    if (xus[3*l]*(wknb+j)->x+xus[1+3*l]*(wknb+j)->y+
                                        xus[2+3*l]*(wknb+j)->z > (wknb+j)->dot)
                                    {
                                        last = j;
                                        break;
                                    }
                                }
                                if (j >= nnei)
                                {
                                    i_ac++;
                                    wkdot[l] = 1;
                                }
                            } /* end of cycle j */
                        }     /* end of cycle l */
                    }
                    else
                    {
                        i_ac  = n_dot;
                        for (l = 0; l < n_dot; l++)
                        {
                            wkdot[l] = 1;
                        }
                    }

                    if (debug)
                    {
                        fprintf(debug, "i_ac=%d, dotarea=%8.3f, aisq=%8.3f\n",
                                i_ac, dotarea, aisq);
                    }

                    a    = aisq*dotarea* (real) i_ac;
                    area = area + a;
                    if (mode & FLAG_ATOM_AREA)
                    {
                        range_check(wkatm[iat], 0, nat);
                        atom_area[wkatm[iat]] = a;
                    }
                    if (mode & FLAG_DOTS)
                    {
                        for (l = 0; l < n_dot; l++)
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
                        dx = 0.; dy = 0.; dz = 0.;
                        for (l = 0; l < n_dot; l++)
                        {
                            if (wkdot[l])
                            {
                                dx = dx+xus[3*l];
                                dy = dy+xus[1+3*l];
                                dz = dz+xus[2+3*l];
                            }
                        }
                        vol = vol+aisq*(dx*(xi-xs)+dy*(yi-ys)+dz*(zi-zs)+ai* (real) i_ac);
                    }

                } /* end of cycle "iat" */
            }     /* end of cycle "ix" */
        }         /* end of cycle "iy" */
    }             /* end of cycle "iz" */

    sfree(wkatm);
    sfree(wkat1);
    sfree(wkdot);
    sfree(wkbox);
    sfree(wknb);
    if (box)
    {
        sfree(x);
    }
    if (mode & FLAG_VOLUME)
    {
        vol           = vol*FOURPI/(3.* (real) n_dot);
        *value_of_vol = vol;
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

    return 0;
}


#if TEST_NSC > 0
#define NAT 2
main () {

    int    i, j, ndots;
    real   co[3*NAT], ra[NAT], area, volume, a, b, c;
    real * dots;
    real * at_area;
    FILE  *fp;


    a  = 1.; c = 0.1;
    fp = fopen("nsc.txt", "w+");
    for (i = 1; i <= NAT; i++)
    {
        j         = i-1;
        co[3*i-3] = j*1*c;
        co[3*i-2] = j*1*c;
        co[3*i-1] = j*1*c;
        /*
           co[3*i-3] = i*1.4;
           co[3*i-2] = 0.;
           co[3*i-1] = 0.;
         */
        /*
           co[3*i-2] = a*0.3;
           a = -a; b=0;
           if (i%3 == 0) b=0.5;
           co[3*i-1] = b;
           ra[i-1] = 2.0;
         */
        ra[i-1] = (1.+j*0.5)*c;
    }
    /*
       if (NSC(co, ra, NAT, 42, NULL, &area,
     */
    if (NSC(co, ra, NAT, 42, NULL, &area,
            NULL, NULL, NULL, NULL))
    {
        ERROR("error in NSC");
    }
    fprintf(fp, "\n");
    fprintf(fp, "area     : %8.3f\n", area);
    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fprintf(fp, "next call\n");
    fprintf(fp, "\n");
    fprintf(fp, "\n");

    if (NSC(co, ra, NAT, 42, FLAG_VOLUME | FLAG_ATOM_AREA | FLAG_DOTS, &area,
            &at_area, &volume,
            &dots, &ndots))
    {
        ERROR("error in NSC");
    }

    fprintf(fp, "\n");
    fprintf(fp, "area     : %8.3f\n", area);
    printf("area     : %8.3f\n", area);
    fprintf(fp, "volume   : %8.3f\n", volume);
    printf("volume   : %8.3f\n", volume);
    fprintf(fp, "ndots    : %8d\n", ndots);
    printf("ndots    : %8d\n", ndots);
    fprintf(fp, "\n");
    for (i = 1; i <= NAT; i++)
    {
        fprintf(fp, "%4d ATOM %7.2f %7.2f %7.2f  ra=%4.1f  area=%8.3f\n",
                i, co[3*i-3], co[3*i-2], co[3*i-1], ra[i-1], at_area[i-1]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "DOTS : %8d\n", ndots);
    for (i = 1; i <= ndots; i++)
    {
        fprintf(fp, "%4d DOTS %8.2f %8.2f %8.2f\n",
                i, dots[3*i-3], dots[3*i-2], dots[3*i-1]);
    }
}
#endif
