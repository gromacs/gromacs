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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _gpp_bondatomtype_h
#define _gpp_bondatomtype_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "macros.h"

typedef struct gpp_bondatomtype *t_bond_atomtype;

extern int get_bond_atomtype_type(char *str,t_bond_atomtype at);
/* Return atomtype corresponding to case-insensitive str
   or NOTSET if not found */

extern char *get_bond_atomtype_name(int nt, t_bond_atomtype at);
/* Return name corresponding to atomtype nt, or NULL if not found */

extern t_bond_atomtype init_bond_atomtype(void);
/* Return a new atomtype structure */

extern void done_bond_atomtype(t_bond_atomtype *at);
/* Free the memory in the structure */

extern void add_bond_atomtype(t_bond_atomtype at,t_symtab *tab,
			      char *name);
/* Add a complete new atom type to an existing atomtype structure */

#endif	/* _gpp_bond_atomtype_h */








