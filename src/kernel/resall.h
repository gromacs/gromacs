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
#include "pdb2gmx.h"
#include "grompp.h"

extern int comprb(const void *a,const void *b);
/* Compare resbond shit */

extern int comprang(const void *a,const void *b);
/* Compare resangle shit */

extern t_resbond *search_rb(char *key,int nrb,t_resbond rb[]);
/* Search in the bond database */

extern t_resang *search_rang(char *key,int nra,t_resang ra[]);
/* Search in the angle database */

extern int comprtp(const void *a,const void *b);
/* Compare routine for sorting and searching for t_restp */

extern t_restp *search_rtp(char *key,int nrtp,t_restp rtp[]);
/* Search for an entry in the rtp database */

extern int icomp(const void *a,const void *b);
/* Compare routine for sorting and searching for t_idihres */

extern t_idihres *search_idih(char *key,int nrdh,t_idihres ires[]);
/* Search for a residue in the improper database */

extern t_atomtype *read_atype(char *adb,t_symtab *tab);

extern int read_resall(char       *resdb,
		       t_restp    **rtp,
		       t_resbond  **rb,
		       t_resang   **ra,
		       t_idihres  **ires,
		       t_atomtype *atype,
		       t_symtab   *tab);

extern void print_resall(FILE *out,
			 int nrtp,
			 t_restp rtp[],
			 t_resbond rb[],
			 t_resang ra[],
			 t_idihres ires[],
			 t_atomtype *atype);

#endif	/* _resall_h */
