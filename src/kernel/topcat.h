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

#ifndef _topcat_h
#define _topcat_h

static char *SRCID_topcat_h = "$Id$";
#ifdef HAVE_IDENT
#ident	"@(#) topcat.h 1.23 9/30/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"

extern void topcat(t_molinfo *dest,int nsrc,t_molinfo src[],
		   int ntab,int *tab,int Nsim,t_simsystem Sims[],
		   bool bEnsemble);
/* If ntab > 0, then molecules will be shuffled over nodes
 * according to tab. If bEnsemble then distance restraints will
 * be added together for ensemble averiging.
 */

extern void mi2top(t_topology *dest,t_molinfo *src);

extern int *mk_shuffle_tab(int nmol,t_molinfo mol[],int nnodes,int *ntab,
			   int Nsim,t_simsystem Sims[],bool bVerbose);
/* Make an array tab (return value) of length *ntab
 * which holds the molecule types
 * which must consecutively be added to the topology
 */
#endif	/* _topcat_h */
