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
 * S  C  A  M  O  R  G
 */

#ifndef	_transfer_h
#define	_transfer_h

#ifdef HAVE_IDENT
#ident	"@(#) transfer.h 1.6 11/23/92"
#endif /* HAVE_IDENT */

extern void linkio_write(int linkno,int nbytes);
extern void linkio_read(int linkno,int nbytes);
extern void linkio_put(int linkno,void *buf,int bufsize);
extern void linkio_get(int linkno,void *buf,int bufsize);

#endif	/* _transfer_h */
