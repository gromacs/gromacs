/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
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
 * GROtesk MACabre and Sinister
 */

#ifndef	_replace_h
#define	_replace_h

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
