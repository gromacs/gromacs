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

#ifndef _utils_h
#define _utils_h

static char *SRCID_utils_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) utils.h 1.11 2/2/97"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"

extern void print_rvec(FILE *log,char *title,rvec vect);

extern void print_rvecs(FILE *log,char *title,int n,rvec vec[]);

#endif	/* _utils_h */
