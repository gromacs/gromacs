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

#ifndef	_shift_h
#define	_shift_h

#ifdef HAVE_IDENT
#ident	"@(#) shift.h 1.6 2/2/97"
#endif /* HAVE_IDENT */
#ifdef HAVE_IDENT
#endif /* HAVE_IDENT */

extern real *mk_shift_tab(int n,real r1,real rc,real dr,real *sfac);
/* Return a table of length n, containing the parabolic
 * shift function from HJC Berendsen
 */


#endif	/* _shift_h */
