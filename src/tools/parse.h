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
 * GROwing Monsters And Cloning Shrimps
 */

#ifndef _parse_h
#define _parse_h

static char *SRCID_parse_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) parse.h 1.17 9/30/97"
#endif /* HAVE_IDENT */
#include <statutil.h>
#include "g_hbond.h"

extern void parse_args(int argc,char *argv[]);

extern void init_topology(char *topology);

extern void user_input(int *nr_groups);

#endif	/* _parse_h */
