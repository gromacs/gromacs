/*
 *       @(#) copyrite.h 1.27 10/15/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.51
 * 
 * Copyright (c) 1990-1996,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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

#ifndef	_copyright_h
#define	_copyright_h

#ifdef HAVE_IDENT
#ident	"@(#) copyright.h 1.10 11/23/92"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include <stdio.h>
  
#define GromacsVersion() "VERSION 2.0b"

static char *CopyrightText[] = {
  "",
  "Copyright (c) 1991-1997",
  "BIOSON Research Institute, Dept. of Biophysical Chemistry",
  "University of Groningen, The Netherlands",
  ""
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
