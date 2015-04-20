/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include <algorithm>

#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"

#define XXXX    0
#define XXYY    1
#define XXZZ    2
#define YYXX    3
#define YYYY    4
#define YYZZ    5
#define ZZXX    6
#define ZZYY    7
#define ZZZZ    8

static void upd_vir(rvec vir, real dvx, real dvy, real dvz)
{
    vir[XX] -= 0.5*dvx;
    vir[YY] -= 0.5*dvy;
    vir[ZZ] -= 0.5*dvz;
}

void calc_vir(int nxf, rvec x[], rvec f[], tensor vir,
              gmx_bool bScrewPBC, matrix box)
{
    int      i, isx;
    double   dvxx = 0, dvxy = 0, dvxz = 0, dvyx = 0, dvyy = 0, dvyz = 0, dvzx = 0, dvzy = 0, dvzz = 0;

    for (i = 0; (i < nxf); i++)
    {
        dvxx += x[i][XX]*f[i][XX];
        dvxy += x[i][XX]*f[i][YY];
        dvxz += x[i][XX]*f[i][ZZ];

        dvyx += x[i][YY]*f[i][XX];
        dvyy += x[i][YY]*f[i][YY];
        dvyz += x[i][YY]*f[i][ZZ];

        dvzx += x[i][ZZ]*f[i][XX];
        dvzy += x[i][ZZ]*f[i][YY];
        dvzz += x[i][ZZ]*f[i][ZZ];

        if (bScrewPBC)
        {
            isx = IS2X(i);
            /* We should correct all odd x-shifts, but the range of isx is -2 to 2 */
            if (isx == 1 || isx == -1)
            {
                dvyy += box[YY][YY]*f[i][YY];
                dvyz += box[YY][YY]*f[i][ZZ];

                dvzy += box[ZZ][ZZ]*f[i][YY];
                dvzz += box[ZZ][ZZ]*f[i][ZZ];
            }
        }
    }

    upd_vir(vir[XX], dvxx, dvxy, dvxz);
    upd_vir(vir[YY], dvyx, dvyy, dvyz);
    upd_vir(vir[ZZ], dvzx, dvzy, dvzz);
}


static void lo_fcv(int i0, int i1,
                   real x[], real f[], tensor vir,
                   int is[], real box[], gmx_bool bTriclinic)
{
    int      i, i3, tx, ty, tz;
    real     xx, yy, zz;
    real     dvxx = 0, dvxy = 0, dvxz = 0, dvyx = 0, dvyy = 0, dvyz = 0, dvzx = 0, dvzy = 0, dvzz = 0;

    if (bTriclinic)
    {
        for (i = i0; (i < i1); i++)
        {
            i3 = DIM*i;
            tx = is[i3+XX];
            ty = is[i3+YY];
            tz = is[i3+ZZ];

            xx    = x[i3+XX]-tx*box[XXXX]-ty*box[YYXX]-tz*box[ZZXX];
            dvxx += xx*f[i3+XX];
            dvxy += xx*f[i3+YY];
            dvxz += xx*f[i3+ZZ];

            yy    = x[i3+YY]-ty*box[YYYY]-tz*box[ZZYY];
            dvyx += yy*f[i3+XX];
            dvyy += yy*f[i3+YY];
            dvyz += yy*f[i3+ZZ];

            zz    = x[i3+ZZ]-tz*box[ZZZZ];
            dvzx += zz*f[i3+XX];
            dvzy += zz*f[i3+YY];
            dvzz += zz*f[i3+ZZ];
        }
    }
    else
    {
        for (i = i0; (i < i1); i++)
        {
            i3 = DIM*i;
            tx = is[i3+XX];
            ty = is[i3+YY];
            tz = is[i3+ZZ];

            xx    = x[i3+XX]-tx*box[XXXX];
            dvxx += xx*f[i3+XX];
            dvxy += xx*f[i3+YY];
            dvxz += xx*f[i3+ZZ];

            yy    = x[i3+YY]-ty*box[YYYY];
            dvyx += yy*f[i3+XX];
            dvyy += yy*f[i3+YY];
            dvyz += yy*f[i3+ZZ];

            zz    = x[i3+ZZ]-tz*box[ZZZZ];
            dvzx += zz*f[i3+XX];
            dvzy += zz*f[i3+YY];
            dvzz += zz*f[i3+ZZ];
        }
    }

    upd_vir(vir[XX], dvxx, dvxy, dvxz);
    upd_vir(vir[YY], dvyx, dvyy, dvyz);
    upd_vir(vir[ZZ], dvzx, dvzy, dvzz);
}

static void lo_fcv2(int i0, int i1,
                    rvec x[], rvec f[], tensor vir,
                    ivec is[], matrix box, gmx_bool bTriclinic)
{
    int      i, gg, tx, ty, tz;
    real     xx, yy, zz;
    real     dvxx = 0, dvxy = 0, dvxz = 0, dvyx = 0, dvyy = 0, dvyz = 0, dvzx = 0, dvzy = 0, dvzz = 0;

    if (bTriclinic)
    {
        for (i = i0, gg = 0; (i < i1); i++, gg++)
        {
            tx = is[gg][XX];
            ty = is[gg][YY];
            tz = is[gg][ZZ];

            xx    = x[i][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
            dvxx += xx*f[i][XX];
            dvxy += xx*f[i][YY];
            dvxz += xx*f[i][ZZ];

            yy    = x[i][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
            dvyx += yy*f[i][XX];
            dvyy += yy*f[i][YY];
            dvyz += yy*f[i][ZZ];

            zz    = x[i][ZZ]-tz*box[ZZ][ZZ];
            dvzx += zz*f[i][XX];
            dvzy += zz*f[i][YY];
            dvzz += zz*f[i][ZZ];
        }
    }
    else
    {
        for (i = i0, gg = 0; (i < i1); i++, gg++)
        {
            tx = is[gg][XX];
            ty = is[gg][YY];
            tz = is[gg][ZZ];

            xx    = x[i][XX]-tx*box[XX][XX];
            dvxx += xx*f[i][XX];
            dvxy += xx*f[i][YY];
            dvxz += xx*f[i][ZZ];

            yy    = x[i][YY]-ty*box[YY][YY];
            dvyx += yy*f[i][XX];
            dvyy += yy*f[i][YY];
            dvyz += yy*f[i][ZZ];

            zz    = x[i][ZZ]-tz*box[ZZ][ZZ];
            dvzx += zz*f[i][XX];
            dvzy += zz*f[i][YY];
            dvzz += zz*f[i][ZZ];
        }
    }

    upd_vir(vir[XX], dvxx, dvxy, dvxz);
    upd_vir(vir[YY], dvyx, dvyy, dvyz);
    upd_vir(vir[ZZ], dvzx, dvzy, dvzz);
}

void f_calc_vir(int i0, int i1, rvec x[], rvec f[], tensor vir,
                t_graph *g, matrix box)
{
    int start, end;

    if (g && (g->nnodes > 0))
    {
        /* Calculate virial for bonded forces only when they belong to
         * this node.
         */
        start = std::max(i0, g->at_start);
        end   = std::min(i1, g->at_end);
#ifdef SAFE
        lo_fcv2(start, end, x, f, vir, g->ishift, box, TRICLINIC(box));
#else
        lo_fcv(start, end, x[0], f[0], vir, g->ishift[0], box[0], TRICLINIC(box));
#endif

        /* If not all atoms are bonded, calculate their virial contribution
         * anyway, without shifting back their coordinates.
         * Note the nifty pointer arithmetic...
         */
        if (start > i0)
        {
            calc_vir(start-i0, x + i0, f + i0, vir, FALSE, box);
        }
        if (end < i1)
        {
            calc_vir(i1-end, x + end, f + end, vir, FALSE, box);
        }
    }
    else
    {
        calc_vir(i1-i0, x + i0, f + i0, vir, FALSE, box);
    }
}
