/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef _systest_h
#define _systest_h

static char *SRCID_systest_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) systest.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

#ifndef MASTER_ID
#define MASTER_ID	2
#endif
#ifndef BUFSIZE
#define BUFSIZE		25000
#endif
#define MEMMARGIN	10000

/*
 * Led usage:
 */
 
#define PARM_LED		0
#define	ACTIVE_LED		1
#define	HOP_LED			2
#define	ERR_LED			3
#define RANDOM_TEST_LED		0
#define WALKBIT_TEST_LED	1
#define WORD_TEST_LED		2
#define RANDOM_SWITCH		0xffff
#define WALKBIT_SWITCH		0x7ffff
#define WORD_SWITCH		0x7ffff

#endif	/* _systest_h */
