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

#ifndef _pgutil_h
#define _pgutil_h

#include "typedefs.h"

extern atom_id search_atom(const char *type,int start,
			   int natoms,t_atom at[],
			   char ** const * anm,
			   const char *bondtype,gmx_bool bDontQuit);
/* Search an atom in array of pointers to strings, starting from start
 * if type starts with '-' then searches backwards from start.
 * bondtype is only used for printing the error/warning string,
 * when bondtype="check" no error/warning is issued.
 * When bDontQuit=FALSE an fatal error is issued, otherwise a warning.
 */

extern void set_at(t_atom *at,real m,real q,int type,int resind);

#endif	/* _pgutil_h */
