/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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
#include "gromacs/math/3dview.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "physics.h"
#include "pbc.h"
#include "vec.h"

#include "gmx_fatal.h"

#define N 4

void m4_op(mat4 m, rvec x, vec4 v)
{
    int i;

    for (i = 0; (i < N); i++)
    {
        v[i] = m[XX][i]*x[XX]+m[YY][i]*x[YY]+m[ZZ][i]*x[ZZ]+m[WW][i];
    }
}

void unity_m4(mat4 m)
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

void print_m4(FILE *fp, const char *s, mat4 A)
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

void print_v4(FILE *fp, char *s, int dim, real *a)
{
    int j;

    if (fp)
    {
        fprintf(fp, "%s: ", s);
        for (j = 0; j < dim; j++)
        {
            fprintf(fp, "%10.5f", a[j]);
        }
        fprintf(fp, "\n");
    }
}

void mult_matrix(mat4 A, mat4 B, mat4 C)
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

void rotate(int axis, real angle, mat4 A)
{
    unity_m4(A);

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

void translate(real tx, real ty, real tz, mat4 A)
{
    unity_m4(A);
    A[3][XX] = tx;
    A[3][YY] = ty;
    A[3][ZZ] = tz;
}

static void set_scale(t_3dview *view, real sx, real sy)
{
    view->sc_x = sx;
    view->sc_y = sy;
}

void calculate_view(t_3dview *view)
{
#define SMALL 1e-6
    mat4 To, Te, T1, T2, T3, T4, T5, N1, D1, D2, D3, D4, D5;
    real dx, dy, dz, l, r;

    /* eye center */
    dx = view->eye[XX];
    dy = view->eye[YY];
    dz = view->eye[ZZ];
    l  = sqrt(dx*dx+dy*dy+dz*dz);
    r  = sqrt(dx*dx+dy*dy);
#ifdef DEBUG
    print_v4(debug, "eye", N, view->eye);
    printf("del: %10.5f%10.5f%10.5f l: %10.5f, r: %10.5f\n", dx, dy, dz, l, r);
#endif
    if (l < SMALL)
    {
        gmx_fatal(FARGS, "Error: Zero Length Vector - No View Specified");
    }
    translate((real)(-view->origin[XX]),
              (real)(-view->origin[YY]), (real)(-view->origin[ZZ]), To);
    translate((real)(-view->eye[XX]),
              (real)(-view->eye[YY]), (real)(-view->eye[ZZ]), Te);

    unity_m4(T2);
    T2[YY][YY] = 0, T2[YY][ZZ] = -1, T2[ZZ][YY] = 1, T2[ZZ][ZZ] = 0;

    unity_m4(T3);
    if (r > 0)
    {
        T3[XX][XX] = -dy/r, T3[XX][ZZ] = dx/r, T3[ZZ][XX] = -dx/r, T3[ZZ][ZZ] = -dy/r;
    }

    unity_m4(T4);
    T4[YY][YY] = r/l, T4[YY][ZZ] = dz/l, T4[ZZ][YY] = -dz/l, T4[ZZ][ZZ] = r/l;

    unity_m4(T5);
    T5[ZZ][ZZ] = -1;

    unity_m4(N1);
    /* N1[XX][XX]=4,N1[YY][YY]=4; */

    mult_matrix(T1, To, view->Rot);
    mult_matrix(D1, Te, T2);
    mult_matrix(D2, T3, T4);
    mult_matrix(D3, T5, N1);
    mult_matrix(D4, T1, D1);
    mult_matrix(D5, D2, D3);

    mult_matrix(view->proj, D4, D5);

#ifdef DEBUG
    print_m4(debug, "T1", T1);
    print_m4(debug, "T2", T2);
    print_m4(debug, "T3", T3);
    print_m4(debug, "T4", T4);
    print_m4(debug, "T5", T5);
    print_m4(debug, "N1", N1);
    print_m4(debug, "Rot", view->Rot);
    print_m4(debug, "Proj", view->proj);
#endif
}

gmx_bool zoom_3d(t_3dview *view, real fac)
{
    real dr;
    real bm, dr1, dr2;
    int  i;

    dr2 = 0;
    for (i = 0; (i < DIM); i++)
    {
        dr   = view->eye[i];
        dr2 += dr*dr;
    }
    dr1 = sqrt(dr2);
    if (fac < 1)
    {
        bm = max(norm(view->box[XX]), max(norm(view->box[YY]), norm(view->box[ZZ])));
        if (dr1*fac < 1.1*bm) /* Don't come to close */
        {
            return FALSE;
        }
    }

    for (i = 0; (i < DIM); i++)
    {
        view->eye[i] *= fac;
    }
    calculate_view(view);
    return TRUE;
}

void init_rotate_3d(t_3dview *view)
{
    real rot = DEG2RAD*15;
    int  i;

    for (i = 0; (i < DIM); i++)
    {
        rotate(i,        rot, view->RotP[i]);
        rotate(i, (real)(-rot), view->RotM[i]);
#ifdef DEBUG
        print_m4(debug, "RotP", view->RotP[i]);
        print_m4(debug, "RotM", view->RotM[i]);
#endif
    }
}


void rotate_3d(t_3dview *view, int axis, gmx_bool bPositive)
{
    int  i, j;
    mat4 m4;

    if (bPositive)
    {
        mult_matrix(m4, view->Rot, view->RotP[axis]);
    }
    else
    {
        mult_matrix(m4, view->Rot, view->RotM[axis]);
    }
    for (i = 0; (i < N); i++)
    {
        for (j = 0; (j < N); j++)
        {
            view->Rot[i][j] = m4[i][j];
        }
    }

    calculate_view(view);
}

void translate_view(t_3dview *view, int axis, gmx_bool bPositive)
{
#ifdef DEBUG
    printf("Translate called\n");
#endif
    if (bPositive)
    {
        view->origin[axis] += view->box[axis][axis]/8;
    }
    else
    {
        view->origin[axis] -= view->box[axis][axis]/8;
    }
    calculate_view(view);
}

void reset_view(t_3dview *view)
{
    int  i;

#ifdef DEBUG
    printf("Reset view called\n");
#endif
    set_scale(view, 4.0, 4.0);
    clear_rvec(view->eye);
    calc_box_center(view->ecenter, view->box, view->origin);
    view->eye[ZZ] = 3.0*max(view->box[XX][XX], view->box[YY][YY]);
    zoom_3d(view, 1.0);
    view->eye[WW] = view->origin[WW] = 0.0;

    /* Initiate the matrix */
    unity_m4(view->Rot);
    calculate_view(view);

    init_rotate_3d(view);
}

t_3dview *init_view(matrix box)
{
    t_3dview *view;
    int       i, j;

    snew(view, 1);

    /* Copy parameters into variables */
    for (i = 0; (i < DIM); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            view->box[i][j] = box[i][j];
        }
    }

    view->ecenter = ecenterDEF;

    reset_view(view);

    return view;
}
