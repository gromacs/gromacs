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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_addconf_c = "$Id$";

#include "vec.h"
#include "macros.h"
#include "smalloc.h"
#include "addconf.h"
#include "gstat.h"


real mydist2(rvec x,rvec y,matrix box)
{
  rvec dx;
  int d;

  for(d=0;(d<DIM);d++) {
    dx[d]=fabs(x[d]-y[d]); 
    while ( dx[d] >= box[d][d] )
      dx[d] -= box[d][d];

    if (fabs(dx[d]-box[d][d])<=dx[d])
      dx[d]-=box[d][d];
  }  


  return iprod(dx,dx);
}

#define MARGE 0.2
bool in_the_box(rvec x,matrix box)
{
  if (((x[0]>=-1 * MARGE)&&(x[0]<=box[0][0]+MARGE))&&
      ((x[1]>=-1 * MARGE)&&(x[1]<=box[1][1]+MARGE))&&
      ((x[2]>=-1 * MARGE)&&(x[2]<=box[2][2]+MARGE)))
    return TRUE;
  else
    return FALSE;
  
}

bool out_the_box(rvec x,matrix box)
{
  if (((x[0]>=MARGE)&&(x[0]<=box[0][0]-MARGE))&&
      ((x[1]>=MARGE)&&(x[1]<=box[1][1]-MARGE))&&
      ((x[2]>=MARGE)&&(x[2]<=box[2][2]-MARGE)))
    return FALSE;
  else
    return TRUE;  
}

void add_conf(t_atoms *atoms_1,rvec *x_1,rvec *v_1,real *r_1,
	      int NTB,matrix box_1,
	      t_atoms *atoms_2,rvec *x_2,rvec *v_2,real *r_2,
	      bool bVerbose,int maxmol)
{
  int  i,j,m,resnr,nresadd,d,k;
  int  *add;
  real d2,vdw2;
  
  if (atoms_2->nr <= 0) {
    fatal_error(0,"Nothing to add");
  }
  
  if (bVerbose)
    fprintf(stderr,"Calculating Overlap...\n");
    
  /* The add array is filled with:
   * 0 for atoms that should NOT be added,
   * 1 for atoms that should be added
   * 2 for atoms that should be added and that are the start of a residue
   */
  snew(add,atoms_2->nr);
  
  /* First set the beginning of residues */
  add[0]  = 2;
  nresadd = 1;
  for (i=1;i<atoms_2->nr;i++) {
    add[i]=1;
    if (atoms_2->atom[i].resnr != atoms_2->atom[i-1].resnr) {
      add[i]  = 2;
      nresadd ++;
    }
  }
  
  
  /* check solvent with solute */
  for(i=0;(i<atoms_1->nr);i++) {
    for(j=0;(j<atoms_2->nr);j++) {
      d2=dist2((x_1)[i],(x_2)[j],box_1);
      vdw2=sqr((r_1)[i] + (r_2)[j]);
      if (d2 < vdw2) {
	if (add[j] == 2)
	  nresadd--;
	add[j] = 0;
      }
      if (add[j]==0) {	
	resnr=atoms_2->atom[j].resnr;
	while((resnr==atoms_2->atom[j-1].resnr) && (j>0))
	  j--;
	if (j<0)
	  warning("j < 0 in addconf");
	while((resnr==atoms_2->atom[j].resnr) && (j<atoms_2->nr)) {
	  if (add[j] == 2)
	    nresadd--;
	  add[j]=0;
	  j++;
	}
      }
    }
  }
  
  /* check solvent with itself */
  for(i=0;(i<atoms_2->nr);i++) {
    bool bAdd = TRUE;

    
    /* remove atoms that are too far away */
    if ((out_the_box(x_2[i],box_1)==TRUE)&&
	(in_the_box(x_2[i],box_1)==FALSE))
      bAdd = FALSE;
    
    /* check only the atoms that are in the border */
    if ( bAdd == TRUE ) {
      if ((out_the_box(x_2[i],box_1)==TRUE)&&
	  (in_the_box(x_2[i],box_1)==TRUE)) {
	
#define EPS 0.00001
	/* check with other border atoms */
	for(j=0;(j<atoms_2->nr);j++) {
	  if (add[j]!=0) {
	    if ( atoms_2->atom[i].resnr!=atoms_2->atom[j].resnr) {
	      if ((in_the_box(x_2[j],box_1)==TRUE)&&
		  (out_the_box(x_2[j],box_1)==TRUE)) {
		d2=mydist2((x_2)[i],(x_2)[j],box_1);
		vdw2=sqr((r_2)[i] + (r_2)[j]);
		if ((d2 < vdw2)||(d2<EPS)) {
		  bAdd = FALSE;
		}
	      }
	    }    
	  }
	}
      }
    }
    
    
    if ( bAdd == FALSE ) {
      for(k=i;((k<atoms_2->nr)&&
	       (atoms_2->atom[k].resnr==atoms_2->atom[i].resnr));k++) {
	if (add[k] == 2)
	  nresadd--;
	add[k] = 0;
      }
      for(k=i;((k>=0)&&
	       (atoms_2->atom[k].resnr==atoms_2->atom[i].resnr));k--) {
	if (add[k] == 2)
	  nresadd--;
	add[k] = 0;
      }
    }
  }
  
  /* add[j]=0 of all the atoms of a molecule containing one or more atoms 
   * add[i]=0
   */
  if (maxmol == 0)
    maxmol = nresadd;
  else
    maxmol = min(maxmol,nresadd);
  fprintf(stderr,"There are %d molecules to add\n",maxmol);
  
  /*add the selected atoms_2 to atoms_1*/  
  for (i=0; (i<atoms_2->nr); i++)
    if (add[i] != 0) {
      if (add[i] == 2) {
	if (maxmol == 0)
	  break;
	atoms_1->nres++;
	maxmol--;
      }
      atoms_1->nr++;
      atoms_1->atomname[atoms_1->nr-1] = atoms_2->atomname[i];
      copy_rvec((x_2)[i],(x_1)[atoms_1->nr-1]);
      copy_rvec((v_2)[i],(v_1)[atoms_1->nr-1]);
      (r_1)[atoms_1->nr-1]   =(r_2)[i];
      atoms_1->atom[atoms_1->nr-1].resnr=atoms_1->nres-1;
      atoms_1->resname[atoms_1->nres-1]=
	atoms_2->resname[atoms_2->atom[i].resnr];
    }
  sfree(add);
}/*add_conf()*/



