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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef	_comtest_h
#define	_comtest_h

#ifdef HAVE_IDENT
#ident	"@(#) comtest.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

extern int master_comtest(char *title,int left,int right,int nprocs,
                          int nbuf,char *fn,int count,int verbose);

extern int slave_comtest(int left,int right,int nprocs);

#endif	/* _comtest_h */
