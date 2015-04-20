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
#include "gmxpre.h"

#include "3dview.h"

#include <math.h>

#include <algorithm>

#include "gromacs/math/3dtransforms.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static void set_scale(t_3dview *view, real sx, real sy)
{
    view->sc_x = sx;
    view->sc_y = sy;
}

static void calculate_view(t_3dview *view)
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
    gmx_vec4_print(debug, "eye", view->eye);
    printf("del: %10.5f%10.5f%10.5f l: %10.5f, r: %10.5f\n", dx, dy, dz, l, r);
#endif
    if (l < SMALL)
    {
        gmx_fatal(FARGS, "Error: Zero Length Vector - No View Specified");
    }
    gmx_mat4_init_translation(-view->origin[XX], -view->origin[YY], -view->origin[ZZ], To);
    gmx_mat4_init_translation(-view->eye[XX], -view->eye[YY], -view->eye[ZZ], Te);

    gmx_mat4_init_unity(T2);
    T2[YY][YY] = 0, T2[YY][ZZ] = -1, T2[ZZ][YY] = 1, T2[ZZ][ZZ] = 0;

    gmx_mat4_init_unity(T3);
    if (r > 0)
    {
        T3[XX][XX] = -dy/r, T3[XX][ZZ] = dx/r, T3[ZZ][XX] = -dx/r, T3[ZZ][ZZ] = -dy/r;
    }

    gmx_mat4_init_unity(T4);
    T4[YY][YY] = r/l, T4[YY][ZZ] = dz/l, T4[ZZ][YY] = -dz/l, T4[ZZ][ZZ] = r/l;

    gmx_mat4_init_unity(T5);
    T5[ZZ][ZZ] = -1;

    gmx_mat4_init_unity(N1);
    /* N1[XX][XX]=4,N1[YY][YY]=4; */

    gmx_mat4_mmul(T1, To, view->Rot);
    gmx_mat4_mmul(D1, Te, T2);
    gmx_mat4_mmul(D2, T3, T4);
    gmx_mat4_mmul(D3, T5, N1);
    gmx_mat4_mmul(D4, T1, D1);
    gmx_mat4_mmul(D5, D2, D3);

    gmx_mat4_mmul(view->proj, D4, D5);

#ifdef DEBUG
    gmx_mat4_print(debug, "T1", T1);
    gmx_mat4_print(debug, "T2", T2);
    gmx_mat4_print(debug, "T3", T3);
    gmx_mat4_print(debug, "T4", T4);
    gmx_mat4_print(debug, "T5", T5);
    gmx_mat4_print(debug, "N1", N1);
    gmx_mat4_print(debug, "Rot", view->Rot);
    gmx_mat4_print(debug, "Proj", view->proj);
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
        bm = std::max(norm(view->box[XX]), std::max(norm(view->box[YY]), norm(view->box[ZZ])));
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

/* Initiates the state of 3d rotation matrices in the structure */
static void init_rotate_3d(t_3dview *view)
{
    real rot = DEG2RAD*15;
    int  i;

    for (i = 0; (i < DIM); i++)
    {
        gmx_mat4_init_rotation(i,  rot, view->RotP[i]);
        gmx_mat4_init_rotation(i, -rot, view->RotM[i]);
#ifdef DEBUG
        gmx_mat4_print(debug, "RotP", view->RotP[i]);
        gmx_mat4_print(debug, "RotM", view->RotM[i]);
#endif
    }
}


void rotate_3d(t_3dview *view, int axis, gmx_bool bPositive)
{
    mat4 m4;

    if (bPositive)
    {
        gmx_mat4_mmul(m4, view->Rot, view->RotP[axis]);
    }
    else
    {
        gmx_mat4_mmul(m4, view->Rot, view->RotM[axis]);
    }
    gmx_mat4_copy(m4, view->Rot);
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
#ifdef DEBUG
    printf("Reset view called\n");
#endif
    set_scale(view, 4.0, 4.0);
    clear_rvec(view->eye);
    calc_box_center(view->ecenter, view->box, view->origin);
    view->eye[ZZ] = 3.0*std::max(view->box[XX][XX], view->box[YY][YY]);
    zoom_3d(view, 1.0);
    view->eye[WW] = view->origin[WW] = 0.0;

    /* Initiate the matrix */
    gmx_mat4_init_unity(view->Rot);
    calculate_view(view);

    init_rotate_3d(view);
}

t_3dview *init_view(matrix box)
{
    t_3dview *view;

    snew(view, 1);
    copy_mat(box, view->box);
    view->ecenter = ecenterDEF;
    reset_view(view);

    return view;
}
