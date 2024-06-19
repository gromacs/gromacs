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
/*! \libinternal
 * \file
 * \brief Defines wrapper functions for higher-level matrix functions
 *
 * \ingroup module_math
 */
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "gromacs/math/nrjac.h"

#include <cmath>
#include <cstdlib>

#include <filesystem>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

template<typename MatrixType>
static inline void do_rotate(MatrixType a, int i, int j, int k, int l, double tau, double s)
{
    double g, h;
    g       = a[i][j];
    h       = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

template<typename MatrixType>
static int jacobi(MatrixType a, const int n, double d[], MatrixType v)
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
        z[ip]         = 0.0;
    }
    int nrot = 0;
    for (i = 1; i <= 50; i++)
    {
        sm = 0.0;
        for (ip = 0; ip < n - 1; ip++)
        {
            for (iq = ip + 1; iq < n; iq++)
            {
                sm += std::abs(a[ip][iq]);
            }
        }
        if (sm == 0.0)
        {
            sfree(z);
            sfree(b);
            return nrot;
        }
        if (i < 4)
        {
            tresh = 0.2 * sm / (n * n);
        }
        else
        {
            tresh = 0.0;
        }
        for (ip = 0; ip < n - 1; ip++)
        {
            for (iq = ip + 1; iq < n; iq++)
            {
                g = 100.0 * std::abs(a[ip][iq]);
                if (i > 4 && std::abs(d[ip]) + g == std::abs(d[ip])
                    && std::abs(d[iq]) + g == std::abs(d[iq]))
                {
                    a[ip][iq] = 0.0;
                }
                else if (std::abs(a[ip][iq]) > tresh)
                {
                    h = d[iq] - d[ip];
                    if (std::abs(h) + g == std::abs(h))
                    {
                        t = (a[ip][iq]) / h;
                    }
                    else
                    {
                        theta = 0.5 * h / (a[ip][iq]);
                        t     = 1.0 / (std::abs(theta) + std::sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                        {
                            t = -t;
                        }
                    }
                    c   = 1.0 / std::sqrt(1 + t * t);
                    s   = t * c;
                    tau = s / (1.0 + c);
                    h   = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j < ip; j++)
                    {
                        do_rotate(a, j, ip, j, iq, tau, s);
                    }
                    for (j = ip + 1; j < iq; j++)
                    {
                        do_rotate(a, ip, j, j, iq, tau, s);
                    }
                    for (j = iq + 1; j < n; j++)
                    {
                        do_rotate(a, ip, j, iq, j, tau, s);
                    }
                    for (j = 0; j < n; j++)
                    {
                        do_rotate(v, j, ip, j, iq, tau, s);
                    }
                    ++nrot;
                }
            }
        }
        for (ip = 0; ip < n; ip++)
        {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    gmx_fatal(FARGS, "Error: Too many iterations in routine JACOBI\n");

    return nrot;
}

void jacobi(double** a, const int numDimensions, double* eigenvalues, double** eigenvectors, int* numRotations)
{
    int numRot = jacobi(a, numDimensions, eigenvalues, eigenvectors);

    if (numRotations)
    {
        *numRotations = numRot;
    }
}

int jacobi(gmx::ArrayRef<gmx::DVec> a, gmx::ArrayRef<double> eigenvalues, gmx::ArrayRef<gmx::DVec> eigenvectors)
{
    GMX_RELEASE_ASSERT(gmx::ssize(a) == DIM, "Size should be 3");
    GMX_RELEASE_ASSERT(gmx::ssize(eigenvalues) == DIM, "Size should be 3");
    GMX_RELEASE_ASSERT(gmx::ssize(eigenvectors) == DIM, "Size should be 3");

    return jacobi(a, DIM, eigenvalues.data(), eigenvectors);
}

int m_inv_gen(const real* m, int n, real* minv)
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
            md[i][j] = m[i * n + j];
        }
    }

    tol = 0;
    for (i = 0; i < n; i++)
    {
        tol += std::abs(md[i][i]);
    }
    tol = 1e-6 * tol / n;

    jacobi(md, n, eig, v, &nrot);

    nzero = 0;
    for (i = 0; i < n; i++)
    {
        if (std::abs(eig[i]) < tol)
        {
            eig[i] = 0;
            nzero++;
        }
        else
        {
            eig[i] = 1.0 / eig[i];
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            s = 0;
            for (k = 0; k < n; k++)
            {
                s += eig[k] * v[i][k] * v[j][k];
            }
            minv[i * n + j] = s;
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
