/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef	_disre_h
#define	_disre_h

#ifdef HAVE_IDENT
#ident	"@(#) disre.h 1.13 2/2/97"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "sysstuff.h"
#include "typedefs.h"

void init_disres(FILE *log,int nbonds,t_inputrec *ir);
/* Initiate local variables, must be called once, nbonds is the number 
 * of iatoms in the ilist of the idef struct
 */

extern real ta_disres(FILE *log,int nbonds,t_iatom fa[],t_iparams *fp,
		      rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
		      matrix box,real lambda,real *dvdlambda,
		      t_mdatoms *md,int ngrp,real egnb[],real egcoul[]);
/* Do the actual calculation */

extern t_drblock *get_drblock(void);

#ifdef CPLUSPLUS
}
#endif

#endif	/* _disre_h */
