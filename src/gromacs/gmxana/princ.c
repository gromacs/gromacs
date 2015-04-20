/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

#include "princ.h"

#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"

static void m_op(matrix mat, rvec x)
{
    rvec xp;
    int  m;

    for (m = 0; (m < DIM); m++)
    {
        xp[m] = mat[m][XX]*x[XX]+mat[m][YY]*x[YY]+mat[m][ZZ]*x[ZZ];
    }
    fprintf(stderr, "x    %8.3f  %8.3f  %8.3f\n", x[XX], x[YY], x[ZZ]);
    fprintf(stderr, "xp   %8.3f  %8.3f  %8.3f\n", xp[XX], xp[YY], xp[ZZ]);
    fprintf(stderr, "fac  %8.3f  %8.3f  %8.3f\n", xp[XX]/x[XX], xp[YY]/x[YY],
            xp[ZZ]/x[ZZ]);
}

#define NDIM 4

#ifdef DEBUG
static void ptrans(char *s, real **inten, real d[], real e[])
{
    int  m;
    real n, x, y, z;
    for (m = 1; (m < NDIM); m++)
    {
        x = inten[m][1];
        y = inten[m][2];
        z = inten[m][3];
        n = x*x+y*y+z*z;
        fprintf(stderr, "%8s %8.3f %8.3f %8.3f, norm:%8.3f, d:%8.3f, e:%8.3f\n",
                s, x, y, z, sqrt(n), d[m], e[m]);
    }
    fprintf(stderr, "\n");
}

void t_trans(matrix trans, real d[], real **ev)
{
    rvec x;
    int  j;
    for (j = 0; (j < DIM); j++)
    {
        x[XX] = ev[1][j+1];
        x[YY] = ev[2][j+1];
        x[ZZ] = ev[3][j+1];
        m_op(trans, x);
        fprintf(stderr, "d[%d]=%g\n", j, d[j+1]);
    }
}
#endif

void principal_comp(int n, atom_id index[], t_atom atom[], rvec x[],
                    matrix trans, rvec d)
{
    int      i, j, ai, m, nrot;
    real     mm, rx, ry, rz;
    double **inten, dd[NDIM], tvec[NDIM], **ev;
#ifdef DEBUG
    real     e[NDIM];
#endif
    real     temp;

    snew(inten, NDIM);
    snew(ev, NDIM);
    for (i = 0; (i < NDIM); i++)
    {
        snew(inten[i], NDIM);
        snew(ev[i], NDIM);
        dd[i] = 0.0;
#ifdef DEBUG
        e[i] = 0.0;
#endif
    }

    for (i = 0; (i < NDIM); i++)
    {
        for (m = 0; (m < NDIM); m++)
        {
            inten[i][m] = 0;
        }
    }
    for (i = 0; (i < n); i++)
    {
        ai           = index[i];
        mm           = atom[ai].m;
        rx           = x[ai][XX];
        ry           = x[ai][YY];
        rz           = x[ai][ZZ];
        inten[0][0] += mm*(sqr(ry)+sqr(rz));
        inten[1][1] += mm*(sqr(rx)+sqr(rz));
        inten[2][2] += mm*(sqr(rx)+sqr(ry));
        inten[1][0] -= mm*(ry*rx);
        inten[2][0] -= mm*(rx*rz);
        inten[2][1] -= mm*(rz*ry);
    }
    inten[0][1] = inten[1][0];
    inten[0][2] = inten[2][0];
    inten[1][2] = inten[2][1];
#ifdef DEBUG
    ptrans("initial", inten, dd, e);
#endif

    for (i = 0; (i < DIM); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            trans[i][m] = inten[i][m];
        }
    }

    /* Call numerical recipe routines */
    jacobi(inten, 3, dd, ev, &nrot);
#ifdef DEBUG
    ptrans("jacobi", ev, dd, e);
#endif

    /* Sort eigenvalues in ascending order */
#define SWAPPER(i)          \
    if (fabs(dd[i+1]) < fabs(dd[i])) {    \
        temp = dd[i];         \
        for (j = 0; (j < NDIM); j++) { tvec[j] = ev[j][i]; } \
        dd[i] = dd[i+1];          \
        for (j = 0; (j < NDIM); j++) { ev[j][i] = ev[j][i+1]; }        \
        dd[i+1] = temp;           \
        for (j = 0; (j < NDIM); j++) { ev[j][i+1] = tvec[j]; }         \
    }
    SWAPPER(0)
    SWAPPER(1)
    SWAPPER(0)
#ifdef DEBUG
    ptrans("swap", ev, dd, e);
    t_trans(trans, dd, ev);
#endif

    for (i = 0; (i < DIM); i++)
    {
        d[i] = dd[i];
        for (m = 0; (m < DIM); m++)
        {
            trans[i][m] = ev[m][i];
        }
    }

    for (i = 0; (i < NDIM); i++)
    {
        sfree(inten[i]);
        sfree(ev[i]);
    }
    sfree(inten);
    sfree(ev);
}

void rotate_atoms(int gnx, atom_id *index, rvec x[], matrix trans)
{
    real   xt, yt, zt;
    int    i, ii;

    for (i = 0; (i < gnx); i++)
    {
        ii        = index ? index[i] : i;
        xt        = x[ii][XX];
        yt        = x[ii][YY];
        zt        = x[ii][ZZ];
        x[ii][XX] = trans[XX][XX]*xt+trans[XX][YY]*yt+trans[XX][ZZ]*zt;
        x[ii][YY] = trans[YY][XX]*xt+trans[YY][YY]*yt+trans[YY][ZZ]*zt;
        x[ii][ZZ] = trans[ZZ][XX]*xt+trans[ZZ][YY]*yt+trans[ZZ][ZZ]*zt;
    }
}

real calc_xcm(rvec x[], int gnx, atom_id *index, t_atom *atom, rvec xcm,
              gmx_bool bQ)
{
    int  i, ii, m;
    real m0, tm;

    clear_rvec(xcm);
    tm = 0;
    for (i = 0; (i < gnx); i++)
    {
        ii = index ? index[i] : i;
        if (atom)
        {
            if (bQ)
            {
                m0 = fabs(atom[ii].q);
            }
            else
            {
                m0 = atom[ii].m;
            }
        }
        else
        {
            m0 = 1;
        }
        tm += m0;
        for (m = 0; (m < DIM); m++)
        {
            xcm[m] += m0*x[ii][m];
        }
    }
    for (m = 0; (m < DIM); m++)
    {
        xcm[m] /= tm;
    }

    return tm;
}

real sub_xcm(rvec x[], int gnx, atom_id *index, t_atom atom[], rvec xcm,
             gmx_bool bQ)
{
    int  i, ii;
    real tm;

    tm = calc_xcm(x, gnx, index, atom, xcm, bQ);
    for (i = 0; (i < gnx); i++)
    {
        ii = index ? index[i] : i;
        rvec_dec(x[ii], xcm);
    }
    return tm;
}

void add_xcm(rvec x[], int gnx, atom_id *index, rvec xcm)
{
    int  i, ii;

    for (i = 0; (i < gnx); i++)
    {
        ii = index ? index[i] : i;
        rvec_inc(x[ii], xcm);
    }
}

void orient_princ(t_atoms *atoms, int isize, atom_id *index,
                  int natoms, rvec x[], rvec *v, rvec d)
{
    int     i, m;
    rvec    xcm, prcomp;
    matrix  trans;

    calc_xcm(x, isize, index, atoms->atom, xcm, FALSE);
    for (i = 0; i < natoms; i++)
    {
        rvec_dec(x[i], xcm);
    }
    principal_comp(isize, index, atoms->atom, x, trans, prcomp);
    if (d)
    {
        copy_rvec(prcomp, d);
    }

    /* Check whether this trans matrix mirrors the molecule */
    if (det(trans) < 0)
    {
        for (m = 0; (m < DIM); m++)
        {
            trans[ZZ][m] = -trans[ZZ][m];
        }
    }
    rotate_atoms(natoms, NULL, x, trans);
    if (v)
    {
        rotate_atoms(natoms, NULL, v, trans);
    }

    for (i = 0; i < natoms; i++)
    {
        rvec_inc(x[i], xcm);
    }
}
