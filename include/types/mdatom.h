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
 * Gyas ROwers Mature At Cryogenic Speed
 */
#ifndef _mdatom_h
#define _mdatom_h

typedef struct {
  real          tmass;
  int           nr;
  real          *massA,*massB,*massT,*invmass;
  real          *chargeA,*chargeB,*chargeT;
  bool          *bPerturbed;
  int           *resnr;
  int           *typeA,*typeB;
  ushort        *ptype;
  ushort        *cTC,*cENER,*cACC,*cFREEZE,*cXTC;
  ushort        *cU1,*cU2,*cU3;
} t_mdatoms;

#endif
