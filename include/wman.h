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

#ifndef _wman_h
#define _wman_h

static char *SRCID_wman_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) wman.h 1.11 10/2/97"
#endif /* HAVE_IDENT */

#include "readinp.h"

enum { eotNone, eotLaTeX, eotHTML, eotNroff, eotAscii, eotJava, eotHelp,
       eotNR };

extern void write_man(FILE *out,int otype,char *program,
		      int nldesc,char *desc[],
		      int nfile,t_filenm fnm[],
		      int npargs,t_pargs pa[],
		      int nbug,char *bugs[],
		      bool bHidden);

extern char *fileopt(ulong flag);
/* Return a string describing the file type in flag.
 * flag should the flag field of a filenm struct.
 */

extern char *check_tty(char *s);
extern char *check_tex(char *s);
extern char *check_html(char *s,char *program);
/* Check LaTeX or HTML strings for codes, and remove them 
 * the program variable may be NULL
 */

#endif	/* _wman_h */


