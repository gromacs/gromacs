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
  return (( (x[0]>=-MARGE) && (x[0]<=box[0][0]+MARGE) ) &&
	  ( (x[1]>=-MARGE) && (x[1]<=box[1][1]+MARGE) ) &&
	  ( (x[2]>=-MARGE) && (x[2]<=box[2][2]+MARGE) ));
}

bool out_the_box(rvec x,matrix box)
{
  return !(( (x[0]>=MARGE) && (x[0]<=box[0][0]-MARGE) )&&
	   ( (x[1]>=MARGE) && (x[1]<=box[1][1]-MARGE) )&&
	   ( (x[2]>=MARGE) && (x[2]<=box[2][2]-MARGE) ));
}

void add_conf(t_atoms *atoms_1,rvec *x_1,rvec *v_1,real *r_1,
	      int NTB,matrix box_1,
	      t_atoms *atoms_2,rvec *x_2,rvec *v_2,real *r_2,
	      bool bVerbose,int maxmol)
{
  enum{Remove,Add,NewRes};
  int  i,j,m,resnr,nresadd,d,k;
  int  *atom_flag;
  real d2,vdw2;
  bool bAdd;
  
  if (atoms_2->nr <= 0) {
    fatal_error(0,"Nothing to add");
  }
  
  if (bVerbose)
    fprintf(stderr,"Calculating Overlap...\n");
  
  /* The atom_flag array is filled with:
   * Remove for atoms that should NOT be added,
   * Add    for atoms that should be added
   * NewRes for atoms that should be added and that are the start of a residue
   */
  snew(atom_flag,atoms_2->nr);
  
  /* First set the beginning of residues */
  atom_flag[0]  = NewRes;
  nresadd = 1;
  for (i=1;i<atoms_2->nr;i++)
    if (atoms_2->atom[i].resnr != atoms_2->atom[i-1].resnr) {
      atom_flag[i] = NewRes;
      nresadd ++;
    } else
      atom_flag[i] = Add;
  
  /* check solvent with solute */
  for(i=0;(i<atoms_1->nr);i++) {
    if ( (i<10) || !((i+1)%10) || (i>atoms_1->nr-10)) 
      fprintf(stderr,"\r%d out of %d atoms checked",i+1,atoms_1->nr);
    for(j=0;(j<atoms_2->nr);j++) {
      d2=dist2((x_1)[i],(x_2)[j],box_1);
      vdw2=sqr((r_1)[i] + (r_2)[j]);
      if (d2 < vdw2) {
	if (atom_flag[j] == NewRes)
	  nresadd--;
	atom_flag[j] = Remove;
      }
      if (atom_flag[j] == Remove) {	
	resnr=atoms_2->atom[j].resnr;
	while((resnr==atoms_2->atom[j-1].resnr) && (j>0))
	  j--;
	if (j<0)
	  warning("j < 0 in addconf");
	while((resnr==atoms_2->atom[j].resnr) && (j<atoms_2->nr)) {
	  if (atom_flag[j] == NewRes)
	    nresadd--;
	  atom_flag[j]=Remove;
	  j++;
	}
      }
    }
  }
  fprintf(stderr,"\n");

  /* check solvent with itself */
#define EPS 0.00001
  for(i=0;(i<atoms_2->nr);i++) {
    if ( (i<10) || !((i+1)%10) || (i>atoms_2->nr-10)) 
      fprintf(stderr,"\r%d out of %d atoms checked",i+1,atoms_2->nr);
    
    /* remove atoms that are too far away */
    bAdd = in_the_box(x_2[i],box_1);
    
    /* check only the atoms that are in the border */
    if ( bAdd && out_the_box(x_2[i],box_1) )
      /* check with other border atoms */
      for(j=0;(j<atoms_2->nr);j++)
	if ( (atom_flag[j] != Remove) && 
	     (atoms_2->atom[i].resnr!=atoms_2->atom[j].resnr) &&
	     in_the_box(x_2[j],box_1) && out_the_box(x_2[j],box_1) ) {
	  d2=mydist2((x_2)[i],(x_2)[j],box_1);
	  vdw2=sqr((r_2)[i] + (r_2)[j]);
	  if ((d2 < vdw2)||(d2<EPS))
	    bAdd = FALSE;
	}
    
    if ( !bAdd ) {
      for(k=i;((k<atoms_2->nr)&&
	       (atoms_2->atom[k].resnr==atoms_2->atom[i].resnr));k++) {
	if (atom_flag[k] == NewRes)
	  nresadd--;
	atom_flag[k] = Remove;
      }
      for(k=i;((k>=0)&&
	       (atoms_2->atom[k].resnr==atoms_2->atom[i].resnr));k--) {
	if (atom_flag[k] == NewRes)
	  nresadd--;
	atom_flag[k] = Remove;
      }
    }
  }
  fprintf(stderr,"\n");
  
  /* atom_flag[j]=Remove of all the atoms of a molecule containing 
   * one or more atoms atom_flag[i]=0
   */
  if (maxmol == 0)
    maxmol = nresadd;
  else
    maxmol = min(maxmol,nresadd);
  fprintf(stderr,"There are %d molecules to add\n",maxmol);
  
  /* add the selected atoms_2 to atoms_1 */
  for (i=0; (i<atoms_2->nr); i++)
    if (atom_flag[i] != Remove) {
      if (atom_flag[i] == NewRes) {
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
  sfree(atom_flag);
}

void orient_mol(t_atoms *atoms,char *indexnm,rvec x[])
{
  int     isize;
  atom_id *index,*simp;
  char    *grpnames;
  real    totmass;
  int     i,m;
  rvec    xcm,angle;
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
  principal_comp(isize,index,atoms->atom,x,trans,angle);
  
  /* Check whether this trans matrix mirrors the molecule */
  if (det(trans) < 0) {
    fprintf(stderr,"Mirroring rotation matrix in Z direction\n");
    for(m=0; (m<DIM); m++)
      trans[ZZ][m] *= -1;
  }  
  rotate_atoms(atoms->nr,simp,x,trans);
  
  if (debug) {
    pr_rvecs(stderr,0,"Rot Matrix",trans,DIM);
    fprintf(stderr,"Det(trans) = %g\n",det(trans));
    
    /* print principal component data */
    fprintf(stderr,"Norm of principal axes before rotation: "
	    "(%.3f, %.3f, %.3f)\n",angle[XX],angle[YY],angle[ZZ]);
    fprintf(stderr,"Totmass = %g\n",totmass);
    principal_comp(isize,index,atoms->atom,x,trans,angle);
    rotate_atoms(atoms->nr,simp,x,trans);
    pr_rvecs(stderr,0,"Rot Matrix",trans,DIM);
  }
  sfree(simp);
  sfree(index);
  sfree(grpnames);
}

