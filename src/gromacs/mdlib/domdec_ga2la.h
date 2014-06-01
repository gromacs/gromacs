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

/* TODO This file should eventually form part of the internals of a
   domain decomposition module */

#ifndef GMX_DOMDEC_GA2LA_H
#define GMX_DOMDEC_GA2LA_H

#include "gmx_ga2la.h"

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Clear all the entries in the ga2la list */
void ga2la_clear(gmx_ga2la_t *ga2la);

gmx_ga2la_t *ga2la_init(int nat_tot, int nat_loc);

/* Set the ga2la entry for global atom a_gl to local atom a_loc and cell. */
void ga2la_set(gmx_ga2la_t *ga2la, int a_gl, int a_loc, int cell);

/* Delete the ga2la entry for global atom a_gl */
void ga2la_del(gmx_ga2la_t *ga2la, int a_gl);

/* Change the local atom for present ga2la entry for global atom a_gl */
void ga2la_change_la(gmx_ga2la_t *ga2la, int a_gl, int a_loc);

/* Returns if the global atom a_gl available locally.
 * Sets the local atom and cell,
 * cell can be larger than the number of zones,
 * in which case it indicates that it is more than one cell away
 * in zone cell - #zones.
 */
gmx_bool ga2la_get(const gmx_ga2la_t *ga2la, int a_gl, int *a_loc, int *cell);

#ifdef __cplusplus
}
#endif

#endif
