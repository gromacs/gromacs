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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */

#ifndef _pgutil_h
#define _pgutil_h

static char *SRCID_pgutil_h = "$Id$";
#ifdef HAVE_IDENT
#ident	"@(#) pgutil.h 1.14 9/30/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"

extern atom_id search_atom(char *type,int start,
			   int natoms,t_atom at[],char **anm[]);
/* Search an atom in array of pointers to strings, starting from start
 * if type starts with '-' then searches backwards from start.
 */

extern void set_at(t_atom *at,real m,real q,int type,int resnr);

#endif	/* _pgutil_h */
