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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _main_h
#define _main_h

static char *SRCID_main_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) main.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "network.h"

extern FILE *stdlog;

extern void open_log(char *fn,t_commrec *cr);
/* Open the log file, if necessary (nprocs > 1) the logfile name is
 * communicated around the ring.
 */

extern t_commrec *init_par(int nprocs,char *argv[]);
/* Initiate the parallel computer. Return the communication record
 * (see network.h). As a side effect the stdlog file is opened.
 */

#endif	/* _main_h */
