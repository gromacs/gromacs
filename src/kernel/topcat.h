/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good gRace! Old Maple Actually Chews Slate
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
/* If ntab > 0, then molecules will be shuffled over processors 
 * according to tab. If bEnsemble then distance restraints will
 * be added together for ensemble averiging.
 */

extern void mi2top(t_topology *dest,t_molinfo *src);

extern int *mk_shuffle_tab(int nmol,t_molinfo mol[],int nprocs,int *ntab,
			   int Nsim,t_simsystem Sims[],bool bVerbose);
/* Make an array tab (return value) of length *ntab
 * which holds the molecule types
 * which must consecutively be added to the topology
 */
#endif	/* _topcat_h */
