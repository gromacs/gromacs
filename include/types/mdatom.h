/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
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
  unsigned short        *ptype;
  unsigned short        *cTC,*cENER,*cACC,*cFREEZE,*cXTC,*cVCM;
  unsigned short        *cU1,*cU2;
} t_mdatoms;

#endif
