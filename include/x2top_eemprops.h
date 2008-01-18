/*
 * $Id$
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

#ifndef _x2top_eemprops_h
#define _x2top_eemprops_h
	
#include <stdio.h>
#include "grompp.h"
	
/* File argument maybe NULL, in which case the default eemprops.dat
 * is opened from the library. If eemtype != -1 eemprops with eemtype
 * equal to eemtype will be read.
 */
extern void *read_eemprops(char *fn,int eemtype);

extern void write_eemprops(FILE *fp,void *eem);

extern int name2eemtype(char *name);

extern int eem_get_numprops(void *eem,int eemtype);

extern int eem_get_index(void *eem,char *resname,char *aname,int eemtype);

extern int eem_get_elem_index(void *eem,int elem,int eemtype);

extern real lo_get_j00(void *eem,int index,real *wj,real qH);

extern real eem_get_j00(void *eem,char *resname,char *aname,
			real *wj,real qH,int eemtype);

extern int eem_get_elem(void *eem,int index);

extern real eem_get_w(void *eem,int index);

extern real eem_get_chi0(void *eem,int index);

extern void eem_set_props(void *eem,int index,real J0,real radius,real chi0);

#endif
