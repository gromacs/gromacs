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
#include "gmxpre.h"

#include "do_fit.h"

#include <math.h>
#include <stdio.h>

#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

real calc_similar_ind(gmx_bool bRho, int nind, atom_id *index, real mass[],
                      rvec x[], rvec xp[])
{
    int  i, j, d;
    real m, tm, xs, xd, rs, rd;

    tm = 0;
    rs = 0;
    rd = 0;
    for (j = 0; j < nind; j++)
    {
        if (index)
        {
            i = index[j];
        }
        else
        {
            i = j;
        }
        m   = mass[i];
        tm += m;
        for (d = 0; d < DIM; d++)
        {
            xd  = x[i][d] - xp[i][d];
            rd += m * sqr(xd);
            if (bRho)
            {
                xs  = x[i][d] + xp[i][d];
                rs += m * sqr(xs);
            }
        }
    }
    if (bRho)
    {
        return 2*sqrt(rd/rs);
    }
    else
    {
        return sqrt(rd/tm);
    }
}

real rmsdev_ind(int nind, atom_id index[], real mass[], rvec x[], rvec xp[])
{
    return calc_similar_ind(FALSE, nind, index, mass, x, xp);
}

real rmsdev(int natoms, real mass[], rvec x[], rvec xp[])
{
    return calc_similar_ind(FALSE, natoms, NULL, mass, x, xp);
}

real rhodev_ind(int nind, atom_id index[], real mass[], rvec x[], rvec xp[])
{
    return calc_similar_ind(TRUE, nind, index, mass, x, xp);
}

real rhodev(int natoms, real mass[], rvec x[], rvec xp[])
{
    return calc_similar_ind(TRUE, natoms, NULL, mass, x, xp);
}

void calc_fit_R(int ndim, int natoms, real *w_rls, rvec *xp, rvec *x, matrix R)
{
    int      c, r, n, j, m, i, irot, s;
    double **omega, **om;
    double   d[2*DIM], xnr, xpc;
    matrix   vh, vk, u;
    real     mn;
    int      index;
    real     max_d;

    if (ndim != 3 && ndim != 2)
    {
        gmx_fatal(FARGS, "calc_fit_R called with ndim=%d instead of 3 or 2", ndim);
    }

    snew(omega, 2*ndim);
    snew(om, 2*ndim);
    for (i = 0; i < 2*ndim; i++)
    {
        snew(omega[i], 2*ndim);
        snew(om[i], 2*ndim);
    }

    for (i = 0; i < 2*ndim; i++)
    {
        d[i] = 0;
        for (j = 0; j < 2*ndim; j++)
        {
            omega[i][j] = 0;
            om[i][j]    = 0;
        }
    }

    /*calculate the matrix U*/
    clear_mat(u);
    for (n = 0; (n < natoms); n++)
    {
        if ((mn = w_rls[n]) != 0.0)
        {
            for (c = 0; (c < ndim); c++)
            {
                xpc = xp[n][c];
                for (r = 0; (r < ndim); r++)
                {
                    xnr      = x[n][r];
                    u[c][r] += mn*xnr*xpc;
                }
            }
        }
    }

    /*construct omega*/
    /*omega is symmetric -> omega==omega' */
    for (r = 0; r < 2*ndim; r++)
    {
        for (c = 0; c <= r; c++)
        {
            if (r >= ndim && c < ndim)
            {
                omega[r][c] = u[r-ndim][c];
                omega[c][r] = u[r-ndim][c];
            }
            else
            {
                omega[r][c] = 0;
                omega[c][r] = 0;
            }
        }
    }

    /*determine h and k*/
    jacobi(omega, 2*ndim, d, om, &irot);
    /*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
     * int     natoms = number of rows and columns
     * real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
     * real       **v = v[0..n-1][0..n-1] contains the vectors in columns
     * int      *irot = number of jacobi rotations
     */

    if (debug && irot == 0)
    {
        fprintf(debug, "IROT=0\n");
    }

    index = 0; /* For the compiler only */

    /* Copy only the first ndim-1 eigenvectors */
    for (j = 0; j < ndim-1; j++)
    {
        max_d = -1000;
        for (i = 0; i < 2*ndim; i++)
        {
            if (d[i] > max_d)
            {
                max_d = d[i];
                index = i;
            }
        }
        d[index] = -10000;
        for (i = 0; i < ndim; i++)
        {
            vh[j][i] = M_SQRT2*om[i][index];
            vk[j][i] = M_SQRT2*om[i+ndim][index];
        }
    }
    if (ndim == 3)
    {
        /* Calculate the last eigenvector as the outer-product of the first two.
         * This insures that the conformation is not mirrored and
         * prevents problems with completely flat reference structures.
         */
        cprod(vh[0], vh[1], vh[2]);
        cprod(vk[0], vk[1], vk[2]);
    }
    else if (ndim == 2)
    {
        /* Calculate the last eigenvector from the first one */
        vh[1][XX] = -vh[0][YY];
        vh[1][YY] =  vh[0][XX];
        vk[1][XX] = -vk[0][YY];
        vk[1][YY] =  vk[0][XX];
    }

    /* determine R */
    clear_mat(R);
    for (r = 0; r < ndim; r++)
    {
        for (c = 0; c < ndim; c++)
        {
            for (s = 0; s < ndim; s++)
            {
                R[r][c] += vk[s][r]*vh[s][c];
            }
        }
    }
    for (r = ndim; r < DIM; r++)
    {
        R[r][r] = 1;
    }

    for (i = 0; i < 2*ndim; i++)
    {
        sfree(omega[i]);
        sfree(om[i]);
    }
    sfree(omega);
    sfree(om);
}

void do_fit_ndim(int ndim, int natoms, real *w_rls, rvec *xp, rvec *x)
{
    int    i, j, m, r, c;
    matrix R;
    rvec   x_old;

    /* Calculate the rotation matrix R */
    calc_fit_R(ndim, natoms, w_rls, xp, x, R);

    /*rotate X*/
    for (j = 0; j < natoms; j++)
    {
        for (m = 0; m < DIM; m++)
        {
            x_old[m] = x[j][m];
        }
        for (r = 0; r < DIM; r++)
        {
            x[j][r] = 0;
            for (c = 0; c < DIM; c++)
            {
                x[j][r] += R[r][c]*x_old[c];
            }
        }
    }
}

void do_fit(int natoms, real *w_rls, rvec *xp, rvec *x)
{
    do_fit_ndim(3, natoms, w_rls, xp, x);
}

void reset_x_ndim(int ndim, int ncm, const atom_id *ind_cm,
                  int nreset, const atom_id *ind_reset,
                  rvec x[], const real mass[])
{
    int  i, m, ai;
    rvec xcm;
    real tm, mm;

    if (ndim > DIM)
    {
        gmx_incons("More than 3 dimensions not supported.");
    }
    tm = 0.0;
    clear_rvec(xcm);
    if (ind_cm != NULL)
    {
        for (i = 0; i < ncm; i++)
        {
            ai = ind_cm[i];
            mm = mass[ai];
            for (m = 0; m < ndim; m++)
            {
                xcm[m] += mm*x[ai][m];
            }
            tm += mm;
        }
    }
    else
    {
        for (i = 0; i < ncm; i++)
        {
            mm = mass[i];
            for (m = 0; m < ndim; m++)
            {
                xcm[m] += mm*x[i][m];
            }
            tm += mm;
        }
    }
    for (m = 0; m < ndim; m++)
    {
        xcm[m] /= tm;
    }

    if (ind_reset != NULL)
    {
        for (i = 0; i < nreset; i++)
        {
            rvec_dec(x[ind_reset[i]], xcm);
        }
    }
    else
    {
        for (i = 0; i < nreset; i++)
        {
            rvec_dec(x[i], xcm);
        }
    }
}

void reset_x(int ncm, const atom_id *ind_cm,
             int nreset, const atom_id *ind_reset,
             rvec x[], const real mass[])
{
    reset_x_ndim(3, ncm, ind_cm, nreset, ind_reset, x, mass);
}
