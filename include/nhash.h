/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _nhash_h
#define _nhash_h

static char *SRCID_nhash_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) nhash.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

#define MAX_LJQQ 997

typedef struct {
  float c6,c12,qq;
} t_ljqq;

extern t_ljqq LJQQ[MAX_LJQQ];

int h_enter(FILE *log,float c6, float c12, float qq);
/* Enters the constants denoted above in the hash-table.
 * Returns the index in LJQQ.
 */

void h_stat(FILE *log);
/* Print statistics for hashing */

#endif	/* _nhash_h */
