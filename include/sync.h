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
 * GROup of MAchos and Cynical Suckers
 */

#ifndef _sync_h
#define _sync_h

static char *SRCID_sync_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) sync.h 1.9 11/23/92"
#endif /* HAVE_IDENT */

#define	SYNC_AVAILABLE	1234
#define	SYNC_DONE	4321

extern void sync_open();
     /*
      * Initialises the synchronisation primimtives by creating the
      * needed semaphores. 
      */
extern void sync_close();
     /*
      * Removes the uses semaphores from the system.
      */
extern void sync_available();
     /*
      * Wait until data is available.
      */
extern void sync_done();
     /*
      * Signal data is handled.
      */
     
#endif	/* _sync_h */
