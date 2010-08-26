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

char *search_rtp(const char *key,int nrtp,t_restp rtp[]);
/* Search for an entry in the rtp database, returns the rtp residue name.
 * A mismatch of one character is allowed, if there is only one nearly
 * matching entry in the database, a warning will be generated.
 */

t_restp *get_restp(const char *rtpname,int nrtp,t_restp rtp[]);
/* Return the entry in the rtp database with rtp name rtpname.
 * Generates a fatal error when rtpname is not found.
 */

gpp_atomtype_t read_atype(const char *ffdir,t_symtab *tab);
/* read atom type database(s) */

void read_resall(char *resdb, int *nrtp,t_restp **rtp, 
		 gpp_atomtype_t atype, t_symtab *tab,
		 gmx_bool bAllowOverrideRTP);
/* read rtp database, append to the existing database */

void print_resall(FILE *out, int nrtp, t_restp rtp[], 
			 gpp_atomtype_t atype);
/* write rtp database */
#ifdef __cplusplus
}
#endif

#endif	/* _resall_h */
