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

#ifndef _ewald_util_h
#define _ewald_util_h

static char *SRCID_ewald_util_h = "$Id$";

#include <math.h>
#include "typedefs.h"
#include "complex.h"


extern real ewald_LRcorrection(FILE *fp,t_nsborder *nsb,
			       t_commrec *cr,t_forcerec *fr,
			       real charge[],t_block *excl,rvec x[],
			       matrix box,rvec mu_tot, real qsum,
			       real surface_eps,matrix lrvir);
/* Calculate the Long range correction to ewald, due to 
 * 1-4 interactions, surface dipole term and charge terms
 */

extern real calc_ewaldcoeff(real rc,real dtol);
/* Determines the Ewald parameter, both for Ewald and PME */
#endif
