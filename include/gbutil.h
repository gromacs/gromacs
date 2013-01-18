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
#include "visibility.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

GMX_LIBGMX_EXPORT
void rotate_conf(int natom, rvec *x, rvec *v, real alfa, real beta, real gamma);
/*rotate() rotates a configuration alfa degrees around the x_axis and beta degrees around the y_axis, *v can be NULL */

void orient(int natom, rvec *x, rvec *v, rvec angle, matrix box);
/*orient() rotates a configuration until the largest atom-atom distance is
 * placed along the z-axis and the second largest distance is placed along
 * the y-axis. Finally the third longest distance is placed along the x-axis
 */

GMX_LIBGMX_EXPORT
void genconf(t_atoms *atoms, rvec *x, rvec *v, real *r, matrix box, ivec n_box);
/*genconf() generates a new configuration by adding boxes*/
GMX_LIBGMX_EXPORT
void gen_box(int NTB, int natoms, rvec *x, matrix box, rvec box_space,
             gmx_bool bCenter);
/* gen_box() generates a box around a configuration, box_space is optional
 * extra space around it. If NTB = 1 then a truncated octahedon will be
 * generated (don't!) if bCenter then coordinates will be centered in the
 * genereated box
 */

#ifdef __cplusplus
}
#endif
