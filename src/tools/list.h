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

#ifndef _list_h
#define _list_h

static char *SRCID_list_h = "$Id$";

#include "hbond.h"
#include "dah.h"
#include <typedefs.h>

class List
{
public:
  List();
  ~List();
private:
  /* list of hydrogen bonds */
  Hbond **hbond;
  
  /* list of all the time frames */
  real  *t;

  /* statistics */
  int  *stat;
  int  *total;

  /* number of hydrogen bonds graph */
  int  **number;

  /* angle and distance distributions */
  real **distance;
  real **angle;
#ifdef ENERGY
  real **energy;
#endif
  
  /* protein hydrogen bond helix */
  int  **helix;

  /* frequency files */
  FILE *fp_freq_intra;
  FILE *fp_freq_inter;
  FILE *fp_freq_all;
  
  /* matrix  */
  char **matr;

 private:
  void assign_type();
  void helical(int i);

  void hydrogenbond_ndx();

  void checktype();
public:
  void print(FILE *output);

  int search(Hbond **dah,int& nr_dah);
  int count(int *nr_dd,atom_id **dd,atom_id **hh,int *nr_aa,atom_id **aa);
  void nosearch(Hbond **dah,int& nr_dah);

  void analyse_init();

  int analyse();

  void dump(t_atoms *atoms);
  int nrtime();
};

#endif


