/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */

#ifndef _toputil_h
#define _toputil_h

static char *SRCID_toputil_h = "$Id$";
#ifdef HAVE_IDENT
#ident	"@(#) toputil.h 1.25 19 Nov 1995"
#endif /* HAVE_IDENT */

#include "grompp.h"

/* UTILITIES */

extern int at2type(char *str, t_atomtype *at);

extern char *type2nm(int nt, t_atomtype *at);

extern void pr_alloc (int extra, t_params *pr);

extern void set_p_string(t_param *p,char *s);

/* INITIATE */

extern void init_plist(t_params plist[]);

extern void init_atomtype (t_atomtype *at);

extern void init_molinfo(t_molinfo *mol);

extern void init_top  (t_topology *top);

extern void done_top(t_topology *top);

/* FREE */
extern void done_block(t_block *block);

extern void done_top(t_topology *top);

extern void done_atom (t_atoms *at);

extern void done_mi(t_molinfo *mi);

/* PRINTING */

extern void print_block(FILE *out,char *szName,char *szIndex, 
			char *szA,t_block *block);

extern void print_atoms(FILE *out,t_atomtype *atype,t_atoms *at,int *cgnr);

extern void print_bondeds(FILE *out,int natoms,directive d,
			  int ftype,int fsubtype,t_params plist[]);

extern void print_excl(FILE *out, int natoms, t_excls excls[]);

#endif	/* _toputil_h */
