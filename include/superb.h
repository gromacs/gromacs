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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _superb_h
#define _superb_h

static char *SRCID_superb_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) superb.h 1.7 2/2/97"
#endif /* HAVE_IDENT */
#include <sysstuff.h>
#include <typedefs.h>

typedef struct {
  t_block *grps;	/* The group members			*/
  char    **name;	/* The group names			*/
  int     *nrdf;	/* Nr of degrees of freedom in a group	*/
  real    *temp;	/* Coupling temperature	per group	*/
  rvec    *acc;		/* Acceleration per group		*/
  tensor  *ekin;	/* Array of energy tensors...		*/
  rvec	  *u;           /* Mean velocities of home particles    */
  atom_id *invgrp;      /* Group number for each atom           */
} t_superblock;

extern t_superblock *init_grps(FILE *log,int left,int right,int pid,int nprocs,
			       char *gfile,bool bMaster);
/* Read a superblock structure from gfile. Do communication if
 * necessary.
 */

#endif	/* _superb_h */
