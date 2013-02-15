/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "macros.h"
#include "vec.h"
#include "gmx_fatal.h"
#include "gstat.h"
#include "pbc.h"
#include "gbutil.h"

static real dist2(t_pbc *pbc, rvec x, rvec y)
{
    rvec dx;

    pbc_dx(pbc, x, y, dx);

    return norm2(dx);
}

real distance_to_z(rvec x)
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

void orient(int natom, rvec *x, rvec *v, rvec angle, matrix box)
{
    real  longest, rij, rzi;
    int   i, j, m, max_i = 0, max_j = 0;
    rvec  origin;
    int   temp;
    real  alfa = 0, beta = 0, gamma = 0;
    t_pbc pbc;

    set_pbc(&pbc, -1, box);

    /*first i am going to look for the longest atom-atom distance*/
    longest = dist2(&pbc, x[0], x[1]);
    i       = 0;
    j       = 1;
    for (i = 0; (i < natom); i++)
    {
        for (j = 0; (j < natom); j++)
        {
            rij = dist2(&pbc, x[i], x[j]);
            if (rij > longest)
            {
                max_i   = i;
                max_j   = j;
                longest = rij;
            }
        }
    }
    /* first check if x[max_i]<x[max_j] else swap*/
    if (x[max_i][2] > x[max_j][2])
    {
        temp  = max_i;
        max_i = max_j;
        max_j = temp;
    }

    /*set the origin to x[i]*/
    for (m = 0; (m < DIM); m++)
    {
        origin[m] = x[max_i][m];
    }
    for (i = 0; (i < natom); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            x[i][m] -= origin[m];
        }
    }

    /* calculate the rotation angles alfa(x_axis) and beta(y_axis)
     * the rotation angles must be calculated clockwise looking
     * along the rotation axis to the origin*
     * alfa (x-axis)
     */
    alfa = atan(x[max_j][ZZ]/x[max_j][YY])-M_PI_2;
    beta = M_PI_2-atan(x[max_j][ZZ]/x[max_j][XX]);
    rotate_conf(natom, x, v, alfa, beta, gamma);

    /* now search the longest distance for rotation along the z_axis */
    longest = distance_to_z(x[0]);
    max_i   = 0;
    for (i = 1; (i < natom); i++)
    {
        rzi = distance_to_z(x[i]);
        if (rzi > longest)
        {
            longest = rzi;
            max_i   = i;
        }
    }
    gamma = atan(x[max_i][YY]/x[max_i][XX])-M_PI_2;
    rotate_conf(natom, x, v, 0, 0, gamma);
    angle[0] = alfa;
    angle[1] = beta;
    angle[2] = gamma;
} /*orient()*/


void genconf(t_atoms *atoms, rvec *x, rvec *v, real *r, matrix box, ivec n_box)
{
    int     i, ix, iy, iz, m, j, imol, offset;
    rvec    delta;
    int     nmol;

    nmol = n_box[XX]*n_box[YY]*n_box[ZZ];

    /*print message*/
    fprintf(stderr, "Generating configuration\n");
    imol = 0;
    for (ix = 0; (ix < n_box[XX]); ix++)
    {
        delta[XX] = ix*box[XX][XX];
        for (iy = 0; (iy < n_box[YY]); iy++)
        {
            delta[YY] = iy*box[YY][YY];
            for (iz = 0; (iz < n_box[ZZ]); iz++)
            {
                delta[ZZ] = iz*box[ZZ][ZZ];
                offset    = imol*atoms->nr;
                for (i = 0; (i < atoms->nr); i++)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        x[offset+i][m] = delta[m]+x[i][m];
                    }
                    if (v)
                    {
                        for (m = 0; (m < DIM); m++)
                        {
                            v[offset+i][m] = v[i][m];
                        }
                    }
                    r[offset+i] = r[i];
                }
                imol++;
            }
        }
    }
    for (i = 1; (i < nmol); i++)
    {
        int offs    = i*atoms->nr;
        int offsres = i*atoms->nres;
        for (j = 0; (j < atoms->nr); j++)
        {
            atoms->atomname[offs+j]                    = atoms->atomname[j];
            atoms->atom[offs+j].resind                 = atoms->atom[j].resind + offsres;
            atoms->resinfo[atoms->atom[offs+j].resind] =
                atoms->resinfo[atoms->atom[j].resind];
            atoms->resinfo[atoms->atom[offs+j].resind].nr += offsres;
        }
    }
    atoms->nr   *= nmol;
    atoms->nres *= nmol;
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            box[j][i] *= n_box[j];
        }
    }
} /*genconf()*/

/*gen_box() generates a box around a configuration*/
void gen_box(int NTB, int natoms, rvec *x, matrix box, rvec box_space,
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

    /*calculate the box size if NTB=1 (truncated octahedron)*/
    if (NTB == 1)
    {
        max_box = box[0][0];
        for (m = 0; (m < DIM); m++)
        {
            max_box = max(max_box, box[m][m]);
        }
        for (m = 0; (m < DIM); m++)
        {
            box[m][m] = max_box;
        }
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


#ifdef DEBUG
    /* print data to check this */
    print_stat(x, natoms, box);
#endif
} /*gen_box()*/
