/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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

#ifndef _steep_h
#define _steep_h

static char *SRCID_steep_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) steep.h 1.6 03 Mar 1996"
#endif /* HAVE_IDENT */
extern real f_max(FILE *log,
		  int left,int right,int nprocs,
		  int start,int end,rvec grad[]);
/* Globally calculate max force */

extern real f_norm(FILE *log,
		  int left,int right,int nprocs,
		  int start,int end,rvec grad[]);
/* Calculates norm of forcee */

extern time_t do_steep(FILE *log,int nfile,t_filenm fnm[],
		       t_parm *parm,t_topology *top,
		       t_groups *grps,t_nsborder *nsb,
		       rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms,
		       tensor ekin,real ener[],
		       t_nrnb nrnb[],
		       bool bVerbose,t_commrec *cr,
		       t_graph *graph);
/* Do steepest descents EM or something like that! */

#endif	/* _steep_h */
