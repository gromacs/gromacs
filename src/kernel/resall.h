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
		       t_atomtype *atype, t_symtab *tab);
/* read rtp database */

extern void print_resall(FILE *out, int bts[], int nrtp, t_restp rtp[], 
			 t_atomtype *atype);
/* write rtp database */

#endif	/* _resall_h */
