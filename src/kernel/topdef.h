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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _topdef_h
#define _topdef_h

static char *SRCID_topdef_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) topdef.h 1.9 11/23/92"
#endif /* HAVE_IDENT */

/* These are the indices for combination rule selection 	*/
#define COMB_GROMOS 	1	/* Gromos rules			*/
#define COMB_EPSSIG 	2	/* Epsilon and Sigma		*/

#endif	/* _topdef_h */
