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

#ifndef _toputil_h
#define _toputil_h

#include "grompp.h"
#include "gpp_atomtype.h"

/* UTILITIES */

extern int name2index(char *str, char ***typenames, int ntypes);

extern void pr_alloc (int extra, t_params *pr);

extern void set_p_string(t_param *p,char *s);

/* INITIATE */

extern void init_plist(t_params plist[]);

extern void init_molinfo(t_molinfo *mol);

extern void init_top  (t_topology *top);

extern void done_top(t_topology *top);

/* FREE */
extern void done_block(t_block *block);

extern void done_top(t_topology *top);

extern void done_atom (t_atoms *at);

extern void done_mi(t_molinfo *mi);

/* PRINTING */

extern void print_blocka(FILE *out,char *szName,char *szIndex, 
			 char *szA,t_blocka *block);

extern void print_atoms(FILE *out,t_atomtype atype,t_atoms *at,int *cgnr);

extern void print_bondeds(FILE *out,int natoms,directive d,
			  int ftype,int fsubtype,t_params plist[]);

extern void print_excl(FILE *out, int natoms, t_excls excls[]);

#endif	/* _toputil_h */
