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
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _replace_h
#define _replace_h

static char *SRCID_replace_h = "$Id$";
#ifdef HAVE_IDENT
#ident	"@(#) replace.h 1.16 10/14/97"
#endif /* HAVE_IDENT */
extern char *replace(char *string,char *search,char *replace);
/* Replace all occurences of 
 * string 'search' in string 'string' by 'replace' 
 */

extern char *replaceww(char *string,char *search,char *replace);
/* Replace all occurences of string 'search' delimited by non-alphanum
 * characters (i.e. whole words) in string 'string' by 'replace' 
 */
#endif	/* _replace_h */
