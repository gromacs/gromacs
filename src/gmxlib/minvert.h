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
 * GROwing Monsters And Cloning Shrimps
 */

#ifndef _minvert_h
#define _minvert_h

static char *SRCID_minvert_h = "$Id$";

#include "typedefs.h"

/* A bunch of routines that works on matrices that run from 1 thru n
 * although they are allocated from 0 thru n
 */
	
extern void mat_mult(int n,real **A,real **B,real **C);

extern real **mk_mat(int n);

extern real **mk_mat2(int nrow,int ncol);

extern void cp_mat(int n,real **src,real **dest);

extern void print_mat(FILE *fp,char *title,int n,real **a,int *indx);
/* index may be NULL */

extern void invert_mat(int n,real **A,real **Ainv);

#endif
