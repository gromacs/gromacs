/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _nmol_h
#define _nmol_h

#include "x11.h"
#include "xutil.h"

extern t_molwin *init_mw(t_x11 *x11,Window Parent,
			 int x,int y,int width,int height,
			 unsigned long fg,unsigned long bg,
			 int ePBC,matrix box);
/* Create the molecule window using the x,y etc. */

extern void map_mw(t_x11 *x11,t_molwin *mw);

extern void z_fill(t_manager *man, real *zz);
extern void create_visibility(t_manager *man);
extern int  compare_obj(const void *a,const void *b);
extern int  filter_vis(t_manager *man);
extern void set_sizes(t_manager *man,real sx,real sy);

extern gmx_bool toggle_hydrogen(t_x11 *x11,t_molwin *mw);
/* Toggle the state of the hydrogen drawing,
 * return the current state
 */

extern void set_bond_type(t_x11 *x11,t_molwin *mw,int bt);
/* Set the state of the atoms drawing. */

extern void set_box_type (t_x11 *x11,t_molwin *mw,int bt);
/* Set the type of box or none (bt = 0)
 */

extern void done_mw(t_x11 *x11,t_molwin *mw);

#endif	/* _nmol_h */
