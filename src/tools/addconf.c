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
 * GRowing Old MAkes el Chrono Sweat
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

static real box_margin;

real max_dist(rvec *x, real *r, int start, int end)
{
  real maxd;
  int i,j;
  
  maxd=0;
  for(i=start; i<end; i++)
    for(j=i+1; j<end; j++)
      maxd=max(maxd,sqrt(distance2(x[i],x[j]))+0.5*(r[i]+r[j]));
  
  return 0.5*maxd;
}

void set_margin(t_atoms *atoms, rvec *x, real *r)
{
  int i,d,start;
/*   char *resname; */

  box_margin=0;
  
  start=0;
  for(i=0; i < atoms->nr; i++) {
    if ( (i+1 == atoms->nr) || 
	 (atoms->atom[i+1].resnr != atoms->atom[i].resnr) ) {
      d=max_dist(x,r,start,i+1);
      if (debug && d>box_margin)
	fprintf(debug,"getting margin from %s: %g\n",
		*(atoms->resname[atoms->atom[i].resnr]),box_margin);
      box_margin=max(box_margin,d);
      start=i+1;
    }
  }
}


bool in_box_plus_margin(rvec x,matrix box)
{
  return ( x[XX]>=-box_margin && x[XX]<=box[XX][XX]+box_margin &&
	   x[YY]>=-box_margin && x[YY]<=box[YY][YY]+box_margin &&
	   x[ZZ]>=-box_margin && x[ZZ]<=box[ZZ][ZZ]+box_margin );
}

bool outside_box_minus_margin(rvec x,matrix box)
{
  return ( x[XX]<box_margin || x[XX]>box[XX][XX]-box_margin ||
	   x[YY]<box_margin || x[YY]>box[YY][YY]-box_margin ||
	   x[ZZ]<box_margin || x[ZZ]>box[ZZ][ZZ]-box_margin );
}

int mark_remove_res(int at, bool *remove, int natoms, t_atom *atom)
{
  int resnr;
  
  resnr = atom[at].resnr;
  while( (at > 0) && (resnr==atom[at-1].resnr) )
    at--;
  while( (at < natoms) && (resnr==atom[at].resnr) ) {
    remove[at]=TRUE;
    at++;
  }
  
  return at;
}

void add_conf(t_atoms *atoms, rvec **x, real **r,  bool bSrenew,
	      int NTB, matrix box,
	      t_atoms *atoms_solvt, rvec *x_solvt, real *r_solvt, 
	      bool bVerbose)
{
  int  i,j,m,prev,resnr,nresadd,d,k;
  int  *atom_flag;
  rvec dx;
  bool *remove;

  if (atoms_solvt->nr <= 0) {
    fprintf(stderr,"WARNING: Nothing to add\n");
    return;
  }
  
  if (bVerbose)
    fprintf(stderr,"Calculating Overlap...\n");
  
  snew(remove,atoms_solvt->nr);
  init_pbc(box,FALSE);
  
  /* set margin around box edges to largest solvent dimension */
  set_margin(atoms_solvt,x_solvt,r_solvt);
  
  /* remove atoms that are far outside the box */
  for(i=0; i<atoms_solvt->nr ;i++)
    if ( !in_box_plus_margin(x_solvt[i],box) )
      i=mark_remove_res(i,remove,atoms_solvt->nr,atoms_solvt->atom);
  
  /* check solvent with solute */
  for(i=0; i<atoms->nr; i++) {
    if ( bVerbose && ( (i<10) || ((i+1)%10) || (i>atoms->nr-10) ) )
      fprintf(stderr,"\r%d out of %d atoms checked",i+1,atoms->nr);
    for(j=0; j<atoms_solvt->nr; j++) 
      /* only check solvent that hasn't been marked for removal already */
      if (!remove[j]) {
	pbc_dx((*x)[i],x_solvt[j],dx);
	if ( norm2(dx) < sqr((*r)[i] + r_solvt[j]) )
	  j=mark_remove_res(j,remove,atoms_solvt->nr,atoms_solvt->atom);
      }
  }
  if (bVerbose)
    fprintf(stderr,"\n");
  
  /* check solvent with itself */
  for(i=0; i<atoms_solvt->nr ;i++) {
    if ( bVerbose && ( (i<10) || !((i+1)%10) || (i>atoms_solvt->nr-10) ) )
      fprintf(stderr,"\rchecking atom %d out of %d",i+1,atoms_solvt->nr);
    
    /* check only atoms that haven't been marked for removal already and
       are close to the box edges */
    if ( !remove[i] && outside_box_minus_margin(x_solvt[i],box) )
      /* check with other border atoms, 
	 only if not already removed and two different molecules (resnr) */
      for(j=i+1; (j<atoms_solvt->nr) && !remove[i]; j++)
	if ( !remove[j] && 
	     atoms_solvt->atom[i].resnr != atoms_solvt->atom[j].resnr &&
	     outside_box_minus_margin(x_solvt[j],box) ) {
	  pbc_dx(x_solvt[i],x_solvt[j],dx);
	  if ( norm2(dx) < sqr(r_solvt[i] + r_solvt[j]) )
	    j=mark_remove_res(j,remove,atoms_solvt->nr,atoms_solvt->atom);
	}
  }
  if (bVerbose)
    fprintf(stderr,"\n");

  /* count how many atoms will be added and make space */
  j=0;
  for (i=0; (i<atoms_solvt->nr); i++)
    if (!remove[i])
      j++;
  if (bSrenew) {
    srenew(atoms->resname,  atoms->nres+atoms_solvt->nres);
    srenew(atoms->atomname, atoms->nr+j);
    srenew(atoms->atom,     atoms->nr+j);
    srenew(*x,              atoms->nr+j);
    srenew(*r,              atoms->nr+j);
  }
  
  /* add the selected atoms_solvt to atoms */
  prev=NOTSET;
  nresadd=0;
  for (i=0; (i<atoms_solvt->nr); i++)
    if (!remove[i]) {
      if ( (prev==NOTSET) || 
	   (atoms_solvt->atom[i].resnr != atoms_solvt->atom[prev].resnr) ) {
	nresadd ++;
	atoms->nres++;
      }
      atoms->nr++;
      atoms->atomname[atoms->nr-1] = atoms_solvt->atomname[i];
      copy_rvec(x_solvt[i],(*x)[atoms->nr-1]);
      (*r)[atoms->nr-1]   = r_solvt[i];
      atoms->atom[atoms->nr-1].resnr = atoms->nres-1;
      atoms->resname[atoms->nres-1] =
	atoms_solvt->resname[atoms_solvt->atom[i].resnr];
      prev=i;
    }
  if (bSrenew)
    srenew(atoms->resname,  atoms->nres+nresadd);
  
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
  fprintf(stderr,"\nSelect group for orientation of molecule:\n");
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

