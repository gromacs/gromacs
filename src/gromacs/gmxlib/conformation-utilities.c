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

#include "conformation-utilities.h"

#include <math.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"

static real dist2(t_pbc *pbc, rvec x, rvec y)
{
    rvec dx;

    pbc_dx(pbc, x, y, dx);

    return norm2(dx);
}

static real distance_to_z(rvec x)
{
    return (sqr(x[XX])+sqr(x[YY]));
} /*distance_to_z()*/

static void low_rotate_conf(int natom, rvec *x, real alfa, real beta, real gamma)
{
    int  i;
    rvec x_old;

    for (i = 0; i < natom; i++)
    {
        copy_rvec(x[i], x_old);
        /*calculate new x[i] by rotation alfa around the x-axis*/
        x[i][XX] =   x_old[XX];
        x[i][YY] =             cos(alfa)*x_old[YY] - sin(alfa)*x_old[ZZ];
        x[i][ZZ] =             sin(alfa)*x_old[YY] + cos(alfa)*x_old[ZZ];
        copy_rvec(x[i], x_old);
        /*calculate new x[i] by rotation beta around the y-axis*/
        x[i][XX] =   cos(beta)*x_old[XX]           + sin(beta)*x_old[ZZ];
        x[i][YY] =                       x_old[YY];
        x[i][ZZ] = -sin(beta)*x_old[XX]           + cos(beta)*x_old[ZZ];
        copy_rvec(x[i], x_old);
        /*calculate new x[i] by rotation gamma around the z-axis*/
        x[i][XX] = x_old[XX]*cos(gamma) - x_old[YY]*sin(gamma);
        x[i][YY] = x_old[XX]*sin(gamma) + x_old[YY]*cos(gamma);
        x[i][ZZ] =                                             x_old[ZZ];
    }
}

static void low_rotate_conf_indexed(int nindex, atom_id *index, rvec *x, real alfa, real beta, real gamma)
{
    int  i;
    rvec x_old;

    for (i = 0; i < nindex; i++)
    {
        copy_rvec(x[index[i]], x_old);
        /*calculate new x[index[i]] by rotation alfa around the x-axis*/
        x[index[i]][XX] =   x_old[XX];
        x[index[i]][YY] =             cos(alfa)*x_old[YY] - sin(alfa)*x_old[ZZ];
        x[index[i]][ZZ] =             sin(alfa)*x_old[YY] + cos(alfa)*x_old[ZZ];
        copy_rvec(x[index[i]], x_old);
        /*calculate new x[index[i]] by rotation beta around the y-axis*/
        x[index[i]][XX] =   cos(beta)*x_old[XX]           + sin(beta)*x_old[ZZ];
        x[index[i]][YY] =                       x_old[YY];
        x[index[i]][ZZ] = -sin(beta)*x_old[XX]           + cos(beta)*x_old[ZZ];
        copy_rvec(x[index[i]], x_old);
        /*calculate new x[index[i]] by rotation gamma around the z-axis*/
        x[index[i]][XX] = x_old[XX]*cos(gamma) - x_old[YY]*sin(gamma);
        x[index[i]][YY] = x_old[XX]*sin(gamma) + x_old[YY]*cos(gamma);
        x[index[i]][ZZ] =                                             x_old[ZZ];
    }
}

void rotate_conf(int natom, rvec *x, rvec *v, real alfa, real beta, real gamma)
{
    if (x)
    {
        low_rotate_conf(natom, x, alfa, beta, gamma);
    }
    if (v)
    {
        low_rotate_conf(natom, v, alfa, beta, gamma);
    }
}

/* Make a new box around a configuration*/
void make_new_box(int natoms, rvec *x, matrix box, rvec box_space,
                  gmx_bool bCenter)
{
    int  i, m;
    rvec xmin, xmax;
    real max_box;

    /*calculate minimum and maximum x[0..DIM-1]*/
    for (m = 0; (m < DIM); m++)
    {
        xmin[m] = xmax[m] = x[0][m];
    }
    for (i = 1; (i < natoms); i++)
    {
        for (m = 0; m < DIM; m++)
        {
            xmin[m] = min(xmin[m], x[i][m]);
            xmax[m] = max(xmax[m], x[i][m]);
        }
    }

    /*calculate the new box sizes for cubic and octahedral ...*/
    for (m = 0; (m < DIM); m++)
    {
        box[m][m] = xmax[m]-xmin[m]+2*box_space[m];
    }

    /*move the molecule to the center of the box*/
    if (bCenter)
    {
        for (i = 0; (i < natoms); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                x[i][m] += 0.5*(box[m][m]-xmin[m]-xmax[m]);
            }
        }
    }
}
