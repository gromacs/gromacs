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
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _vcm_h
#define _vcm_h

static char *SRCID_vcm_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) vcm.h 1.9 9/29/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "typedefs.h"
	
extern void calc_vcm(FILE *log,int homenr,int start,
		     real mass[],rvec v[],rvec vcm);

extern void do_stopcm(FILE *log,int homenr,int start,
		      rvec v[],rvec mvcm,real tm);

extern void check_cm(FILE *log,rvec mvcm,real tm);

/* remove global rotation of system by fitting to structure of nstcomm
   steps ago */
extern void do_stoprot(FILE *log, int natoms, rvec box, rvec x[], 
		       real mass[]);

#endif /* _vcm_h */
