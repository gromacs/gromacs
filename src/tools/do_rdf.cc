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
static char *SRCID_do_rdf_cc = "$Id$";

#include <math.h>
#include "maths.h"
#include "vec.h"
#include "smalloc.h"
#include "macros.h"
#include "txtdump.h"
#include "pbc.h"
#include "statutil.h"
#include "string2.h"
#include "statutil.h"
#include "do_rdf.h"
	
#define RESOLUTION 0.005   /* 0.005 nm resolution */

c_rdf::c_rdf(int n1,char *gnm1,int n2,char *gnm2,bool bCM)
{
  nx1=n1;
  nm1=strdup(gnm1);
  nx2=n2;
  nm2=strdup(gnm2);
  bCoM=bCM;
  
  ng=0;
}

c_rdf::~c_rdf()
{
  free(inv_segvol);
}

void c_rdf::init_box(matrix box)
{
  int   i;
  real  boxmin;

  boxmin=min(norm(box[XX]),min(norm(box[YY]),norm(box[ZZ])));
  hbox2=boxmin/2.0;
  if (ng == 0) {
    ng=(int)(hbox2/RESOLUTION)+1;
    init_graph();
  }
  segsize=RESOLUTION;
  rseg_2=1.0/sqr(segsize);
  
  inv_segvol=(real *)calloc(ng,sizeof(inv_segvol[0]));
  
  /* Calculate volume of the spere */
  for(i=0; (i<ng); i++)
    inv_segvol[i]=sphere_vol((i+1)*segsize);
    
  /* Now only count the segment */
  for(i=ng-1; (i>0); i--)
    inv_segvol[i]-=inv_segvol[i-1];

  /* Now take the inverse */
  for(i=0; (i<ng); i++)
    inv_segvol[i]=1.0/inv_segvol[i];
}

void c_rdf::add(real cutoff,t_block *excl,matrix box,rvec x[], 
		atom_id index1[],atom_id index2[])
{
  real    hb2;
  int     i,j,j0,j1;
  int     ix,jx;
  atom_id *index3;
  bool    *bExcl;
  rvec    dx,xcm;
  real    r2,cut2;

  hb2=sqr(hbox2);
  cut2=sqr(cutoff);
  
  /* Must init pbc every step because of pressure coupling */
  init_pbc(box,FALSE);
  j0=0;
  if (index2) {
    j1=nx2;
    index3=index2;
  }
  else {
    j1=nx1;
    index3=index1;
  }
  if (bCoM) {
    clear_rvec(xcm);
    for(i=0; (i < nx1); i++) {
      ix=index1[i];
      rvec_inc(xcm,x[ix]);
    }
    
    for(j=0; (j<DIM); j++)
      xcm[j]/=nx1;
#ifdef DEBUGRDF
    fprintf(stderr,"xcm=(%10.3f, %10.3f, %10.3f)\n",xcm[0],xcm[1],xcm[2]);
#endif      
    for(j=j0; (j < j1); j++) {
      jx=index3[j];
      rvec_sub(x[jx],xcm,dx);
      r2=iprod(dx,dx);
      if ((r2 < hb2) && (r2 > cut2)) 
	insert(dx,r2);
    }
  }
  else {
    for(i=0; (i < nx1); i++) {
      ix=index1[i];
      if (excl) {
	bExcl=(bool *)calloc(excl->nr,sizeof(*bExcl));
	for(j=excl->index[ix]; (j<excl->index[ix+1]); j++)
	  bExcl[excl->a[j]]=TRUE;
      } else
      if (!index2)
	j0=i+1;
      for(j=j0; (j < j1); j++) {
	jx=index3[j];
	if (!excl || !bExcl[jx]) {
	  pbc_dx(box,x[ix],x[jx],dx);
	  r2=iprod(dx,dx);
	  if ((r2 < hb2) && (r2 > cut2)) 
	    insert(dx,r2);
	}
      }
      free(bExcl);
    }
  }
}

void c_rdf::calc(char *fn,real cutoff,t_block *excl,
		 atom_id index1[],atom_id index2[])
{
  real         t;
  int          j,natoms,status;
  rvec         *x0;
  matrix       box;
  
  if ((natoms=read_first_x(&status,fn,&t,&x0,box))==0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }
  init_box(box);

  j=0;
  do {
    add(cutoff,excl,box,x0,index1,index2);
    j++;
  } while (read_next_x(status,&t,natoms,x0,box));
  fprintf(stderr,"\n");
  
  close_trj(status);
  
  sfree(x0);
}

real sphere_vol(real r)
{
  return (4.0*M_PI/3.0)*(r*r*r);
}
