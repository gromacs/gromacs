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

#ifndef _resall_h
#define _resall_h

#include "typedefs.h"
#include "hackblock.h"
#include "gpp_atomtype.h"
#include "grompp.h"

#ifdef __cplusplus
extern "C" {
#endif

extern t_restp *search_rtp(const char *key,int nrtp,t_restp rtp[]);
/* Search for an entry in the rtp database */

extern gpp_atomtype_t read_atype(const char *ffdir,bool bAddCWD,t_symtab *tab);
/* read atom type database(s) */

extern void read_resall(char *resdb, int *nrtp,t_restp **rtp, 
			gpp_atomtype_t atype, t_symtab *tab,
			bool bAllowOverrideRTP);
/* read rtp database, append to the existing database */

extern void print_resall(FILE *out, int nrtp, t_restp rtp[], 
			 gpp_atomtype_t atype);
/* write rtp database */
#ifdef __cplusplus
}
#endif

#endif	/* _resall_h */
