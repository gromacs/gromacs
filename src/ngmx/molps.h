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

#ifndef _molps_h
#define _molps_h

static char *SRCID_molps_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) molps.h 1.10 9/30/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "manager.h"
	
extern void ps_draw_mol(FILE *ps,t_manager *man);
/* Draw molecules to a postscript file */

#endif	/* _molps_h */
