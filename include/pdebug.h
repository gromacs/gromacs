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

#ifndef _pdebug_h
#define _pdebug_h

static char *SRCID_pdebug_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) pdebug.h 1.7 2/2/97"
#endif /* HAVE_IDENT */
#ifdef HAVE_IDENT
#endif /* HAVE_IDENT */
extern void p_debug(char *s,char *file,int line);

#ifdef  DEBUG
#define PDEBUG(s) p_debug(s,__FILE__,__LINE__)
#else
#define PDEBUG(s)
#endif

extern void pma(FILE *log,char *file,int line);

#ifdef DEBUG
#define PMA() pma(log,__FILE__,__LINE__)
#else
#define PMA()
#endif

#endif	/* _pdebug_h */
