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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _add_par_h
#define _add_par_h

#include "typedefs.h"
#include "pdb2top.h"

extern void add_param(t_params *ps, int ai, int aj, real *c, char *s);

extern void add_imp_param(t_params *ps, int ai, int aj, int ak, int al,
			  real c0, real c1, char *s);
			  
extern void add_dih_param(t_params *ps,int ai,int aj,int ak,int al,
			  real c0, real c1, real c2, char *s);

extern void add_cmap_param(t_params *ps,int ai,int aj,int ak,int al,int am,
						   char *s);

extern void add_vsite2_atoms(t_params *ps, int ai, int aj, int ak);

extern void add_vsite3_atoms(t_params *ps, int ai, int aj, int ak, int al, 
			   gmx_bool bSwapParity);

extern void add_vsite2_param(t_params *ps, int ai, int aj, int ak, real c0);

extern void add_vsite3_param(t_params *ps, int ai, int aj, int ak, int al, 
			   real c0, real c1);

extern void add_vsite4_atoms(t_params *ps, int ai, int aj, int ak, int al, 
			   int am);

extern int search_jtype(t_restp *rp,char *name,gmx_bool bFirstRes);

#endif	/* _add_par_h */
