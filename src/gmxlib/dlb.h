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
 * GROwing Monsters And Cloning Shrimps
 */

#ifndef _dlb_h
#define _dlb_h

static char *SRCID_dlb_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) dlb.h 1.3 31 Jan 1995"
#endif /* HAVE_IDENT */
#include "typedefs.h"

extern void count_nb(t_commrec *cr,t_nsborder *nsb,t_block *cgs,int nns,
		     int nlr,t_idef *idef,int ngner);

#endif	/* _dlb_h */
