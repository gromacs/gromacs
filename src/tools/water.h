/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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

#ifndef _water_h
#define _water_h

static char *SRCID_water_h = "$Id$";

#include "g_hbond.h"
#include <typedefs.h>

class Water {
 public:
  Water(Hbond *hbond);
  ~Water();
 private:
  atom_id nr;
  bool    insert;
 public:
  bool exist();
  void print(FILE *fp);
};

void water_list_init(atom_id *w,int nrw);

#endif










