/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef _nmol_h
#define _nmol_h

#include "manager.h"
#include "x11.h"
#include "xutil.h"

extern t_molwin *init_mw(t_x11 *x11, Window Parent,
                         int x, int y, int width, int height,
                         unsigned long fg, unsigned long bg,
                         int ePBC, matrix box);
/* Create the molecule window using the x,y etc. */

extern void map_mw(t_x11 *x11, t_molwin *mw);

extern void z_fill(t_manager *man, real *zz);
extern int  compare_obj(const void *a, const void *b);
extern int  filter_vis(t_manager *man);
extern void set_sizes(t_manager *man);

extern bool toggle_hydrogen(t_x11 *x11, t_molwin *mw);
/* Toggle the state of the hydrogen drawing,
 * return the current state
 */

extern void set_bond_type(t_x11 *x11, t_molwin *mw, int bt);
/* Set the state of the atoms drawing. */

extern void set_box_type (t_x11 *x11, t_molwin *mw, int bt);
/* Set the type of box or none (bt = 0)
 */

extern void done_mw(t_x11 *x11, t_molwin *mw);

#endif  /* _nmol_h */
