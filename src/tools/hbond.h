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

#ifndef _hbond_h
#define _hbond_h

static char *SRCID_hbond_h = "$Id$";

#include <typedefs.h>
#include "g_hbond.h"

class Hbond
{
 public:
  Hbond(atom_id dd,atom_id aa,atom_id hh);
  Hbond(Hbond *that);
  ~Hbond();
 private:
  /* the atom numbers of the donor, acceptor and hydrogen atom of this hbond */
  atom_id d,a,h;

  /* the time in ps of the first and the last occurence of this hbond */
  real first;
  real last;

  /* the total number of times this hbond has been observed */
  int  freq;

  /* distance of this hbond , selected only */
  real *dist;

  /* water distances for insert only */
  atom_id *waid;
  real    *d_dist;
  real    *a_dist;
  real    *h_dist;

 public:
  t_helix_hb hhb;
  t_type     type;

  real distance2();
  real angle();

#ifdef ENERGY
  real energy();
#endif

  bool exist(void);

  bool insert(void);
  void insert_print(int t,FILE *fp);
  void selected_print(int t,FILE *fp);
  
  void print(FILE *output);

  /* set type of this hbond */
  void settype(atom_id left,atom_id right);

  
  /* analyse frequency etc */
  void analyse_init(void);
  void analyse(void);

  int  compare(Hbond *that);

  char *hb_name(t_atoms *atoms);  

  void ndx(t_type t,FILE *fp);

};

void water_list_init(atom_id *wp,int nrwp);

#endif






