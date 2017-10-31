/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2016,2017, by the GROMACS development team, led by
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

#include "3dtransforms.h"

#include <stdio.h>

#include <cmath>

#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"

#define N 4

void gmx_mat4_copy(mat4 a, mat4 b)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            b[i][j] = a[i][j];
        }
    }
}

void gmx_mat4_transform_point(mat4 m, const rvec x, vec4 v)
{
    int i;

    for (i = 0; (i < N); i++)
    {
        v[i] = m[XX][i]*x[XX]+m[YY][i]*x[YY]+m[ZZ][i]*x[ZZ]+m[WW][i];
    }
}

void gmx_mat4_mmul(mat4 A, mat4 B, mat4 C)
{
    int i, j, k;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i][j] = 0;
            for (k = 0; (k < N); k++)
            {
                A[i][j] += B[i][k]*C[k][j];
            }
        }
    }
}

void gmx_mat4_init_unity(mat4 m)
{
    int i, j;

    for (i = 0; (i < N); i++)
    {
        for (j = 0; (j < N); j++)
        {
            if (i == j)
            {
                m[i][j] = 1.0;
            }
            else
            {
                m[i][j] = 0.0;
            }
        }
    }
}

void gmx_mat4_init_rotation(int axis, real angle, mat4 A)
{
    gmx_mat4_init_unity(A);

    switch (axis)
    {
        case XX:
            A[YY][YY] =  cos(angle);
            A[YY][ZZ] = -sin(angle);
            A[ZZ][YY] =  sin(angle);
            A[ZZ][ZZ] =  cos(angle);
            break;
        case YY:
            A[XX][XX] =  cos(angle);
            A[XX][ZZ] =  sin(angle);
            A[ZZ][XX] = -sin(angle);
            A[ZZ][ZZ] =  cos(angle);
            break;
        case ZZ:
            A[XX][XX] =  cos(angle);
            A[XX][YY] = -sin(angle);
            A[YY][XX] =  sin(angle);
            A[YY][YY] =  cos(angle);
            break;
        default:
            gmx_fatal(FARGS, "Error: invalid axis: %d", axis);
    }
}

void gmx_mat4_init_translation(real tx, real ty, real tz, mat4 A)
{
    gmx_mat4_init_unity(A);
    A[3][XX] = tx;
    A[3][YY] = ty;
    A[3][ZZ] = tz;
}

void gmx_mat4_print(FILE *fp, const char *s, mat4 A)
{
    int i, j;

    if (fp)
    {
        fprintf(fp, "%s: ", s);
        for (i = 0; i < N; i++)
        {
            fprintf(fp, "\t");
            for (j = 0; j < N; j++)
            {
                fprintf(fp, "%10.5f", A[i][j]);
            }
            fprintf(fp, "\n");
        }
    }
}

void gmx_vec4_print(FILE *fp, const char *s, vec4 a)
{
    int j;

    if (fp)
    {
        fprintf(fp, "%s: ", s);
        for (j = 0; j < N; j++)
        {
            fprintf(fp, "%10.5f", a[j]);
        }
        fprintf(fp, "\n");
    }
}
