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

#ifndef _main_h
#define _main_h

static char *SRCID_main_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) main.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "network.h"

extern FILE *stdlog;
extern int  gmx_parallel; /* 1 when running in parallel */

extern char *par_fn(char *base,int ftp,t_commrec *cr);
/* Add processor id in the filename right before the extension */

extern void open_log(char *fn,t_commrec *cr);
/* Open the log file, if necessary (nprocs > 1) the logfile name is
 * communicated around the ring.
 */

extern t_commrec *init_par(int *argc,char ***argv_ptr);
/* Initiate the parallel computer. Return the communication record
 * (see network.h). The command line arguments are communicated so that they can be
 * parsed on each processor.
 * Arguments are the number of command line arguments, and a pointer to the
 * array of argument strings.
 */

#endif	/* _main_h */
