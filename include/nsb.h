/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _nsb_h
#define _nsb_h

static char *SRCID_nsb_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern void calc_nsbshift(FILE *fp,t_nsborder *nsb);
/* Calculates the shift and bshift variables */

extern void calc_nsb(FILE *fp,t_block *cgs,int nprocs,
		     t_nsborder *nsb,int nstDlb);
/* Calculate which blocks of charge groups should be calculated,
 * depending on processor number.
 */

extern void print_nsb(FILE *fp,char *title,t_nsborder *nsb);
/* Print the values to file */

/*extern bool cg_avail(int icg,int jcg,t_nsborder *nsb);*/
/* Determine whether a charge group jcg is in the domain of icg:
 * this is necessary for searching on a grid, to avoid double
 * pairs of interactions.
 */

/*extern int  cg_index(int icg,t_nsborder *nsb);*/
/* Perform a modulo calculation giving the correct cg index */

#endif
