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
 * Great Red Owns Many ACres of Sand 
 */
#include <stdio.h>
#include "typedefs.h"
#include "force.h"
#include "lrutil.h"

int main(int argc,char *argv[])
{
  t_forcerec *fr;
  rvec box;
  
  fr=mk_forcerec();
  fr->r1 = 0.6;
  fr->rc = 0.9;
  fr->eeltype = eelTWIN;
  box[XX]=box[YY]=box[ZZ]=1.0;
  
  set_LRconsts(stdout,fr->r1,fr->rc,box,fr);

  make_tables(fr);
  
  return 0;
}
