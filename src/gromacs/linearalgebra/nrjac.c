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

#include "nrjac.h"

#include <math.h>

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static gmx_inline
void do_rotate(double **a, int i, int j, int k, int l, double tau, double s)
{
    double g, h;
    g       = a[i][j];
    h       = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
    int    j, i;
    int    iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

    snew(b, n);
    snew(z, n);
    for (ip = 0; ip < n; ip++)
    {
        for (iq = 0; iq < n; iq++)
        {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }
    for (ip = 0; ip < n; ip++)
    {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    *nrot = 0;
    for (i = 1; i <= 50; i++)
    {
        sm = 0.0;
        for (ip = 0; ip < n-1; ip++)
        {
            for (iq = ip+1; iq < n; iq++)
            {
                sm += fabs(a[ip][iq]);
            }
        }
        if (sm == 0.0)
        {
            sfree(z);
            sfree(b);
            return;
        }
        if (i < 4)
        {
            tresh = 0.2*sm/(n*n);
        }
        else
        {
            tresh = 0.0;
        }
        for (ip = 0; ip < n-1; ip++)
        {
            for (iq = ip+1; iq < n; iq++)
            {
                g = 100.0*fabs(a[ip][iq]);
                if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
                    && fabs(d[iq])+g == fabs(d[iq]))
                {
                    a[ip][iq] = 0.0;
                }
                else if (fabs(a[ip][iq]) > tresh)
                {
                    h = d[iq]-d[ip];
                    if (fabs(h)+g == fabs(h))
                    {
                        t = (a[ip][iq])/h;
                    }
                    else
                    {
                        theta = 0.5*h/(a[ip][iq]);
                        t     = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0)
                        {
                            t = -t;
                        }
                    }
                    c         = 1.0/sqrt(1+t*t);
                    s         = t*c;
                    tau       = s/(1.0+c);
                    h         = t*a[ip][iq];
                    z[ip]    -= h;
                    z[iq]    += h;
                    d[ip]    -= h;
                    d[iq]    += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j < ip; j++)
                    {
                        do_rotate(a, j, ip, j, iq, tau, s);
                    }
                    for (j = ip+1; j < iq; j++)
                    {
                        do_rotate(a, ip, j, j, iq, tau, s);
                    }
                    for (j = iq+1; j < n; j++)
                    {
                        do_rotate(a, ip, j, iq, j, tau, s);
                    }
                    for (j = 0; j < n; j++)
                    {
                        do_rotate(v, j, ip, j, iq, tau, s);
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip = 0; ip < n; ip++)
        {
            b[ip] +=  z[ip];
            d[ip]  =  b[ip];
            z[ip]  =  0.0;
        }
    }
    gmx_fatal(FARGS, "Error: Too many iterations in routine JACOBI\n");
}

int m_inv_gen(real **m, int n, real **minv)
{
    double **md, **v, *eig, tol, s;
    int      nzero, i, j, k, nrot;

    snew(md, n);
    for (i = 0; i < n; i++)
    {
        snew(md[i], n);
    }
    snew(v, n);
    for (i = 0; i < n; i++)
    {
        snew(v[i], n);
    }
    snew(eig, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            md[i][j] = m[i][j];
        }
    }

    tol = 0;
    for (i = 0; i < n; i++)
    {
        tol += fabs(md[i][i]);
    }
    tol = 1e-6*tol/n;

    jacobi(md, n, eig, v, &nrot);

    nzero = 0;
    for (i = 0; i < n; i++)
    {
        if (fabs(eig[i]) < tol)
        {
            eig[i] = 0;
            nzero++;
        }
        else
        {
            eig[i] = 1.0/eig[i];
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            s = 0;
            for (k = 0; k < n; k++)
            {
                s += eig[k]*v[i][k]*v[j][k];
            }
            minv[i][j] = s;
        }
    }

    sfree(eig);
    for (i = 0; i < n; i++)
    {
        sfree(v[i]);
    }
    sfree(v);
    for (i = 0; i < n; i++)
    {
        sfree(md[i]);
    }
    sfree(md);

    return nzero;
}
