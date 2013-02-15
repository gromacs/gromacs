/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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

#ifndef _INTERACTION_CONST_
#define _INTERACTION_CONST_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /* VdW */
    real rvdw;
    real sh_invrc6; /* For shifting the LJ potential */

    /* type of electrostatics (defined in enums.h) */
    int  eeltype;

    /* Coulomb */
    real rcoulomb;

    /* Cut-off */
    real rlist;
    real rlistlong;

    /* PME/Ewald */
    real ewaldcoeff;
    real sh_ewald;   /* For shifting the Ewald potential */

    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r;
    real epsfac;

    /* Constants for reaction-field or plain cut-off */
    real epsilon_rf;
    real k_rf;
    real c_rf;

    /* Force/energy interpolation tables, linear in force, quadratic in V */
    real  tabq_scale;
    int   tabq_size;
    /* Coulomb force table, size of array is tabq_size (when used) */
    real *tabq_coul_F;
    /* Coulomb energy table, size of array is tabq_size (when used) */
    real *tabq_coul_V;
    /* Coulomb force+energy table, size of array is tabq_size*4,
       entry quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
       this is used with single precision x86 SIMD for aligned loads */
    real *tabq_coul_FDV0;
} interaction_const_t;

#ifdef __cplusplus
}
#endif

#endif /* _INTERACTION_CONST_ */
