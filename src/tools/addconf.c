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
#include "princ.h"
#include "rdgroup.h"
#include "txtdump.h"
#include "pbc.h"

#define MARGIN 0.2
bool in_box_plus_margin(rvec x,matrix box)
{
  return ( ( (x[0]>=-MARGIN) && (x[0]<=box[0][0]+MARGIN) ) &&
	   ( (x[1]>=-MARGIN) && (x[1]<=box[1][1]+MARGIN) ) &&
	   ( (x[2]>=-MARGIN) && (x[2]<=box[2][2]+MARGIN) ) );
}

bool outside_box_minus_margin(rvec x,matrix box)
{
  return ( ( (x[0]<MARGIN) && (x[0]>box[0][0]-MARGIN) ) &&
	   ( (x[1]<MARGIN) && (x[1]>box[1][1]-MARGIN) ) &&
	   ( (x[2]<MARGIN) && (x[2]>box[2][2]-MARGIN) ) );
}

void add_conf(t_atoms *atoms_1,rvec *x_1,rvec *v_1,real *r_1,
	      int NTB,matrix box_1,
	      t_atoms *atoms_2,rvec *x_2,rvec *v_2,real *r_2,
	      bool bVerbose)
{
  int  i,j,m,prev,resnr,nresadd,d,k;
  int  *atom_flag;
  rvec dx;
  bool bRemove;
  bool *remove;
  
  if (atoms_2->nr <= 0) {
    fprintf(stderr,"WARNING: Nothing to add\n");
    return;
  }
  
  if (bVerbose)
    fprintf(stderr,"Calculating Overlap...\n");
  
  snew(remove,atoms_2->nr);
  init_pbc(box_1,FALSE);
  
  /* check solvent with solute */
  for(i=0;(i<atoms_1->nr);i++) {
    if ( bVerbose && ( (i<10) || ((i+1)%10) || (i>atoms_1->nr-10) ) )
      fprintf(stderr,"\r%d out of %d atoms checked",i+1,atoms_1->nr);
    for(j=0;(j<atoms_2->nr);j++) 
      /* only check solvent that hasn't been marked for removal already */
      if (!remove[j]) {
	pbc_dx(x_1[i],x_2[j],dx);
	remove[j] = norm2(dx) < sqr(r_1[i] + r_2[j]);
	if (remove[j]) {	
	  resnr=atoms_2->atom[j].resnr;
	  while( (j > 0) && (resnr==atoms_2->atom[j-1].resnr) )
	    j--;
	  for( ; (j < atoms_2->nr) && (resnr==atoms_2->atom[j].resnr); j++)
	    remove[j]=TRUE;
	}
      }
  }
  if (bVerbose)
    fprintf(stderr,"\n");

  /* check solvent with itself */
  for(i=0;(i<atoms_2->nr);i++) {
    if ( bVerbose && ( (i<10) || !((i+1)%10) || (i>atoms_2->nr-10) ) )
      fprintf(stderr,"\rchecking atom %d out of %d",i+1,atoms_2->nr);
    
    /* remove atoms that are too far away */
    remove[i] = !in_box_plus_margin(x_2[i],box_1);
    
    /* check only the atoms that are in the border */
    if ( !remove[i] && outside_box_minus_margin(x_2[i],box_1) )
      /* check with other border atoms, 
	 only if not already removed and two different molecules (resnr) */
      for(j=0; (j<atoms_2->nr) && !remove[i]; j++)
	if ( !remove[j] && 
	     atoms_2->atom[i].resnr != atoms_2->atom[j].resnr &&
	     in_box_plus_margin(x_2[j],box_1) && 
	     outside_box_minus_margin(x_2[j],box_1) ) {
	  pbc_dx(x_2[i],x_2[j],dx);
	  remove[i] = norm2(dx) < sqr(r_2[i] + r_2[j]);
	}
    
    if ( remove[i] ) {
      k=i;
      while( (k>=0) && (atoms_2->atom[k].resnr==atoms_2->atom[i].resnr) )
	k--;
      for(; (k < atoms_2->nr) && 
	    (atoms_2->atom[k].resnr==atoms_2->atom[i].resnr); k++)
	remove[k] = TRUE;
    }
  }
  if (bVerbose)
    fprintf(stderr,"\n");
  
  /* add the selected atoms_2 to atoms_1 */
  prev=NOTSET;
  nresadd=0;
  for (i=0; (i<atoms_2->nr); i++)
    if (!remove[i]) {
      if ( (prev==NOTSET) || 
	   (atoms_2->atom[i].resnr != atoms_2->atom[prev].resnr) ) {
	nresadd ++;
	atoms_1->nres++;
      }
      atoms_1->nr++;
      atoms_1->atomname[atoms_1->nr-1] = atoms_2->atomname[i];
      copy_rvec(x_2[i],x_1[atoms_1->nr-1]);
      copy_rvec(v_2[i],v_1[atoms_1->nr-1]);
      r_1[atoms_1->nr-1]   =r_2[i];
      atoms_1->atom[atoms_1->nr-1].resnr=atoms_1->nres-1;
      atoms_1->resname[atoms_1->nres-1]=
	atoms_2->resname[atoms_2->atom[i].resnr];
      prev=i;
    }
  if (bVerbose)
    fprintf(stderr,"Added %d molecules\n",nresadd);
  
  sfree(remove);
}

void orient_mol(t_atoms *atoms,char *indexnm,rvec x[], rvec *v)
{
  int     isize;
  atom_id *index,*simp;
  char    *grpnames;
  real    totmass;
  int     i,m;
  rvec    xcm,prcomp;
  matrix  trans;
  
  /* Make an index for principal component analysis */
  fprintf(stderr,"Select group for orientation of molecule:\n");
  get_index(atoms,indexnm,1,&isize,&index,&grpnames);
  snew(simp,atoms->nr);
  for(i=0; (i<atoms->nr); i++) {
    simp[i]=i;
    atoms->atom[i].m=1;
  }
  totmass = sub_xcm(x,atoms->nr,simp,atoms->atom,xcm,FALSE);
  principal_comp(isize,index,atoms->atom,x,trans,prcomp);
  
  /* Check whether this trans matrix mirrors the molecule */
  if (det(trans) < 0) {
    if (debug)
      fprintf(stderr,"Mirroring rotation matrix in Z direction\n");
    for(m=0; (m<DIM); m++)
      trans[ZZ][m] = -trans[ZZ][m];
  }  
  rotate_atoms(atoms->nr,simp,x,trans);
  
  if (debug) {
    pr_rvecs(stderr,0,"Rot Matrix",trans,DIM);
    fprintf(stderr,"Det(trans) = %g\n",det(trans));
    
    /* print principal component data */
    fprintf(stderr,"Norm of principal axes before rotation: "
	    "(%.3f, %.3f, %.3f)\n",prcomp[XX],prcomp[YY],prcomp[ZZ]);
    fprintf(stderr,"Totmass = %g\n",totmass);
    principal_comp(isize,index,atoms->atom,x,trans,prcomp);
    rotate_atoms(atoms->nr,simp,x,trans);
    if (v) 
      rotate_atoms(atoms->nr,simp,v,trans);
    pr_rvecs(stderr,0,"Rot Matrix",trans,DIM);
  }
  sfree(simp);
  sfree(index);
  sfree(grpnames);
}

