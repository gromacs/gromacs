/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */

#ifndef _sync_h
#define _sync_h

static char *SRCID_sync_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
