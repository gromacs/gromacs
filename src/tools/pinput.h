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
 * GROup of MAchos and Cynical Suckers
 */
#ifndef _pinput_h
#define _pinput_h

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
