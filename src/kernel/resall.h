/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */

#ifndef _resall_h
#define _resall_h

static char *SRCID_resall_h = "$Id$";
#ifdef HAVE_IDENT
#ident	"@(#) resall.h 1.16 9/30/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "hackblock.h"
#include "grompp.h"

extern t_restp *search_rtp(char *key,int nrtp,t_restp rtp[]);
/* Search for an entry in the rtp database */

extern t_atomtype *read_atype(char *adb,t_symtab *tab);
/* read atom type database */

extern int read_resall(char *resdb, int bts[], t_restp **rtp, 
		       t_atomtype *atype, t_symtab *tab, bool *bAlldih, int *nrexcl);
/* read rtp database */

extern void print_resall(FILE *out, int bts[], int nrtp, t_restp rtp[], 
			 t_atomtype *atype, bool bAlldih, int nrexcl);
/* write rtp database */

#endif	/* _resall_h */
