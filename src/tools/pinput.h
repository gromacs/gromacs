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
 * S  C  A  M  O  R  G
 */

#ifndef _pinput_h
#define _pinput_h

static char *SRCID_pinput_h = "$Id$";

#include "typedefs.h"
#include "string2.h"

enum { ptMC, ptREC, ptPTRJ, ptNR };

typedef struct {
  real step;
  real tol;
  real v0;
  char base[STRLEN];
  char recomb[STRLEN];
  char gamma[STRLEN];
  int  funct;
  int  nsteps;
  int  nframes;
  int  nskip;
  int  nSel;
  int  nev;
} t_pinp;
	
extern void read_inp(char *fnin,char *fnout,t_pinp *p);

#endif
