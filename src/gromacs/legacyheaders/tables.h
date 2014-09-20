/*
 * This file is part of the GROMACS molecular simulation package.
 *
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

#ifndef _tables_h
#define _tables_h
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef double (*real_space_grid_contribution_computer)(double, double);
/* Function pointer used to tell table_spline3_fill_ewald_lr whether it
 * should calculate the grid contribution for electrostatics or LJ.
 */

void table_spline3_fill_ewald_lr(real                                 *table_F,
                                 real                                 *table_V,
                                 real                                 *table_FDV0,
                                 int                                   ntab,
                                 double                                dx,
                                 real                                  beta,
                                 real_space_grid_contribution_computer v_lr);
/* Fill tables of ntab points with spacing dr with the ewald long-range
 * (mesh) force.
 * There are three separate tables with format FDV0, F, and V.
 * This function interpolates the Ewald mesh potential contribution
 * with coefficient beta using a quadratic spline.
 * The force can then be interpolated linearly.
 */

real ewald_spline3_table_scale(const interaction_const_t *ic);
/* Return the scaling for the Ewald quadratic spline tables. */

double v_q_ewald_lr(double beta, double r);
/* Return the real space grid contribution for Ewald*/

double v_lj_ewald_lr(double beta, double r);
/* Return the real space grid contribution for LJ-Ewald*/

#ifdef __cplusplus
}
#endif

#endif  /* _tables\_h */
