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
 * Green Red Orange Magenta Azure Cyan Skyblue
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
