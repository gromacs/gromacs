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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */

#ifndef	_pdebug_h
#define	_pdebug_h

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
