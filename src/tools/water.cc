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
 * GROtesk MACabre and Sinister
 */
static char *SRCID_water_cc = "$Id$";

#include "water.h"
#include <typedefs.h>
#include <pbc.h>
#include <vec.h>

static atom_id *w;
static int     nrw;


real dist(matrix box,rvec v_1,rvec v_2)
{
  rvec dx;
  pbc_dx(box,v_1,v_2,dx);
  return (norm(dx));
}


atom_id find_nearest_atom(atom_id atom,atom_id *water,int nrwater,
			  rvec *x,matrix box)
{
  int i;
  real min_distance=1000;
  atom_id dummy=0;
  for(i=0;(i<nrwater);i++)
    if ((water[i]!=atom)&&(dist(box,x[water[i]],x[atom])<min_distance)) {
      min_distance=dist(box,x[water[i]],x[atom]);
      dummy=water[i];
    }
  return(dummy);
}

void find_atoms(atom_id atom,real rcut,
		atom_id **fw, int& nrfw,
		rvec *x,matrix box)
{
  int i;
  nrfw=0;
  for(i=0;(i<nrw);i++)
    if (dist(box,x[w[i]],x[atom])<rcut) {
      (*fw) = (atom_id *)realloc(*fw,++nrfw*sizeof(atom_id));
      (*fw)[nrfw-1]=w[i];
    }
  if (nrfw==0) {
    (*fw)=(atom_id *)malloc(sizeof(atom_id));
    (*fw)[0]=find_nearest_atom(atom,w,nrw,x,box);
    nrfw=1;
  }
}



Water::Water(Hbond *hbond)
{
  atom_id *water_index=NULL;
  int water_size=0;
  find_atoms(a,rcut,&water_index,water_size,x,box);
  nr=find_nearest_atom(d,water_index,water_size,x,box);
  d_dist=dist(box,x[nr],x[d]);
  a_dist=dist(box,x[nr],x[a]);
  h_dist=dist(box,x[nr],x[h]);
  distan=dist(box,x[d],x[a]);
  if ((a_dist<rcut) && (d_dist<rcut))
    insert=TRUE;
  else
    insert=FALSE;
}

Water::~Water()
{
}

bool Water::exist()
{
  return insert;
}


void water_list_init(atom_id *wp,int nrwp)
{
  w = wp;
  nrw = nrwp;
}

void Water::print(FILE *fp)
{
  fprintf(fp,"%8.3f %8.3f %8.3f %8.3f %5d",distan,d_dist,h_dist,a_dist,nr);
}




