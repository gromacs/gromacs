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
 * Grunge ROck MAChoS
 */

#ifndef _copyrite_h
#define _copyrite_h

static char *SRCID_copyrite_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) copyright.h 1.10 11/23/92"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif
  
#include <stdio.h>
  
/* Has to be a function, so we can get version number from autoconf */   
char *GromacsVersion(void);
  
  
static char *CopyrightText[] = {
  "Copyright (c) 1991-2001, University of Groningen, The Netherlands"
};

static char *GPLText[] = {
  "This program is free software; you can redistribute it and/or",
  "modify it under the terms of the GNU General Public License",
  "as published by the Free Software Foundation; either version 2",
  "of the License, or (at your option) any later version."
};


extern void pr_difftime(FILE *out,double dt);

void CopyRight(FILE *out,char *szProgram);

extern char *bromacs(void);

extern char *cool_quote(void);

extern int be_cool(void);
/* Return TRUE when the user is COOL, FALSE otherwise */

extern void thanx(FILE *fp);

enum { eCITEGMX, eCITEBATH, eCITESHAKE, eCITESETTLE, eCITESOR, 
       eCITEDISRE, eCITERF, eCITELINCS, eCITENR };

extern void please_cite(FILE *fp,char *key);
/* Print a message asking to cite something... */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _copyright_h */
