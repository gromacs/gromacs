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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_pbc_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "vec.h"
#include "maths.h"
#include "main.h"
#include "pbc.h"
#include "smalloc.h"
#include "txtdump.h"
#include "fatal.h"

/*****************************************
 *     PERIODIC BOUNDARY CONDITIONS
 *****************************************/

/* Global variables */
static bool bInit=FALSE;
static bool bTriclinic,bSupported;
static rvec gl_fbox,gl_hbox,gl_mhbox,*tric_vec=NULL;
static matrix gl_box;
static int ntric_vec;
static real sure_dist2;

char *check_box(matrix box)
{
  char *ptr;

  if ((box[XX][YY]!=0) || (box[XX][ZZ]!=0) || (box[YY][ZZ]!=0))
    ptr = "Only triclinic boxes with the first vector parallel to the x-axis and the second vector in the xy-plane are supported.";
  else if ((fabs(box[YY][XX])+fabs(box[ZZ][XX]) > 1.001*box[XX][XX]) ||
	   (fabs(box[ZZ][YY]) > 1.001*box[YY][YY]))
    ptr = "Triclinic box is too skewed.";
  else
    ptr = NULL;
  
  return ptr;
}

void init_pbc(matrix box,bool bTruncOct)
{
  static int nalloc=0;
  int  i,j,k,d;
  real diagonal2;
  rvec try;
  char *ptr;

  ptr = check_box(box);
  if (ptr) {
    fprintf(stderr,"Warning: %s Can not fix pbc.\n",ptr);
    bSupported = FALSE;
  } else {
    bSupported = TRUE;

    for(i=0; (i<DIM); i++) {
      gl_fbox[i]  =  box[i][i];
      gl_hbox[i]  =  gl_fbox[i]*0.5;
      gl_mhbox[i] = -gl_hbox[i];
    }
    bTriclinic = TRICLINIC(box);
    if (bTriclinic) {
      copy_mat(box,gl_box);
      /* When a 'shifted' distance is within this number, it is the shortest
       * possible distance of all shifts.
       */
      sure_dist2 = 0.25*min(norm2(box[XX]),
			    min(norm2(box[YY]),norm2(box[ZZ])));
      /* Make shift vectors, assuming the box is not very skewed */
      diagonal2 = norm2(gl_fbox);
      ntric_vec = 0;
      for(i=-2; i<=2; i++)
	for(j=-2; j<=2; j++)
	  for(k=-2; k<=2; k++)
	    if ((i!=0) || (j!=0) || (k!=0)) {
	      for(d=0; d<DIM; d++)
		try[d] = i*box[XX][d]+j*box[YY][d]+k*box[ZZ][d];
	      if (norm2(try) < diagonal2) {
		if (ntric_vec >= nalloc) {
		  nalloc+=20;
		  srenew(tric_vec,nalloc);
		}
	      copy_rvec(try,tric_vec[ntric_vec]);
	      ntric_vec++;
	      }
	    }
    }
  }
  bInit   = TRUE;
}

void pbc_dx(rvec x1, rvec x2, rvec dx)
{
  int i,j;
  rvec dx_start,try;
  real d2min,d2try;

  if (!bInit)
    fatal_error(0,"pbc_dx called before init_pbc");
  rvec_sub(x1,x2,dx);
  if (bSupported) {
    if (bTriclinic) {
      for(i=DIM-1; i>=0; i--)
	if (dx[i] > gl_hbox[i])
	  for (j=i; j>=0; j--)
	    dx[j] -= gl_box[i][j];
	else if (dx[i] <= gl_mhbox[i])
	  for (j=i; j>=0; j--)
	    dx[j] += gl_box[i][j];
      /* dx is the distance in a rectangular box */
      copy_rvec(dx,dx_start);
      d2min = norm2(dx);
      /* now try all possible shifts */
      i=0;
      while ((d2min > sure_dist2) && (i < ntric_vec)) {
	rvec_add(dx_start,tric_vec[i],try);
	d2try = norm2(try);
	if (d2try < d2min) {
	  copy_rvec(try,dx);
	  d2min = d2try;
	}
	i++;
      }
    } else {
      for(i=0; i<DIM; i++)
	if (dx[i] > gl_hbox[i])
	  dx[i] -= gl_fbox[i];
	else if (dx[i] <= gl_mhbox[i])
	  dx[i] += gl_fbox[i];
    }
  }
}

bool image_rect(ivec xi,ivec xj,ivec box_size,real rlong2,int *shift,real *r2)
{
  int 	m,t;
  int 	dx,b,b_2;
  real  dxr,rij2;

  rij2=0.0;
  t=0;
  for(m=0; (m<DIM); m++) {
    dx=xi[m]-xj[m];
    t*=DIM;
    b=box_size[m];
    b_2=b/2;
    if (dx < -b_2) {
      t+=2;
      dx+=b;
    }
    else if (dx > b_2)
      dx-=b;
    else
      t+=1;
    dxr=dx;
    rij2+=dxr*dxr;
    if (rij2 >= rlong2) 
      return FALSE;
  }
  
  *shift = t;
  *r2 = rij2;
  return TRUE;
}

bool image_cylindric(ivec xi,ivec xj,ivec box_size,real rlong2,
		     int *shift,real *r2)
{
  int 	m,t;
  int 	dx,b,b_2;
  real  dxr,rij2;

  rij2=0.0;
  t=0;
  for(m=0; (m<DIM); m++) {
    dx=xi[m]-xj[m];
    t*=DIM;
    b=box_size[m];
    b_2=b/2;
    if (dx < -b_2) {
      t+=2;
      dx+=b;
    }
    else if (dx > b_2)
      dx-=b;
    else
      t+=1;

    dxr=dx;
    rij2+=dxr*dxr;
    if (m < ZZ) {
      if (rij2 >= rlong2) 
	return FALSE;
    }
  }
  
  *shift = t;
  *r2 = rij2;
  return TRUE;
}

void calc_shifts(matrix box,rvec box_size,rvec shift_vec[],bool bTruncOct)
{
  int k,l,m,d,n,test;
  
  init_pbc(box,bTruncOct);
  for (m=0; (m<DIM); m++)
    box_size[m]=box[m][m];
  
  n=0;
  for(k = -D_BOX; k <= D_BOX; k++) 
    for(l = -D_BOX; l <= D_BOX; l++) 
      for(m = -D_BOX; m <= D_BOX; m++,n++) {
	test=XYZ2IS(k,l,m);
	if (n != test) 
	  fprintf(stdlog,"n=%d, test=%d\n",n,test);
	for(d = 0; d < DIM; d++)
	  shift_vec[n][d]=k*box[XX][d] + l*box[YY][d] + m*box[ZZ][d];
      }
}

void calc_cgcm(FILE *log,int cg0,int cg1,t_block *cgs,
	       rvec pos[],rvec cg_cm[])
{
  int  icg,ai,k,k0,k1,d;
  rvec cg;
  real nrcg,inv_ncg;
  atom_id *cga,*cgindex;
  
#ifdef DEBUG
  fprintf(log,"Calculating centre of geometry for charge groups %d to %d\n",
	  cg0,cg1);
#endif
  cga     = cgs->a;
  cgindex = cgs->index;
  
  /* Compute the center of geometry for all charge groups */
  for(icg=cg0; (icg<cg1); icg++) {
    k0      = cgindex[icg];
    k1      = cgindex[icg+1];
    nrcg    = k1-k0;
    if (nrcg == 1) {
      ai = cga[k0];
      copy_rvec(pos[ai],cg_cm[icg]);
    }
    else {
      inv_ncg = 1.0/nrcg;
      
      clear_rvec(cg);
      for(k=k0; (k<k1); k++)  {
	ai     = cga[k];
	for(d=0; (d<DIM); d++)
	  cg[d] += pos[ai][d];
      }
      for(d=0; (d<DIM); d++)
	cg_cm[icg][d] = inv_ncg*cg[d];
    }
  }
}

void put_charge_groups_in_box(FILE *log,int cg0,int cg1,bool bTruncOct,
			      matrix box,rvec box_size,t_block *cgs,
			      rvec pos[],rvec shift_vec[],rvec cg_cm[])
			      
{
  int  icg,ai,k,k0,k1,d,e;
  rvec cg;
  real nrcg,inv_ncg;
  atom_id *cga,*cgindex;
  bool bTric;

#ifdef DEBUG
  fprintf(log,"Putting cgs %d to %d in box\n",cg0,cg1);
#endif
  cga     = cgs->a;
  cgindex = cgs->index;

  bTric = TRICLINIC(box);

  for(icg=cg0; (icg<cg1); icg++) {
    /* First compute the center of geometry for this charge group */
    k0      = cgindex[icg];
    k1      = cgindex[icg+1];
    nrcg    = k1-k0;
    inv_ncg = 1.0/nrcg;
    
    clear_rvec(cg);
    for(k=k0; (k<k1); k++)  {
      ai     = cga[k];
      for(d=0; d<DIM; d++)
	cg[d] += inv_ncg*pos[ai][d];
    }
    /* Now check pbc for this cg */
    if (bTric) {
      for(d=DIM-1; d>=0; d--) {
	while(cg[d] < 0) {
	  for(e=d; e>=0; e--) {
	    cg[e] += box[d][e];
	    for(k=k0; (k<k1); k++) 
	      pos[cga[k]][e] += box[d][e];
	  }
	}
	while(cg[d] >= box[d][d]) {
	  for(e=d; e>=0; e--) {
	    cg[e] -= box[d][e];
	    for(k=k0; (k<k1); k++) 
	      pos[cga[k]][e] -= box[d][e];
	  }
	}
	cg_cm[icg][d] = cg[d];
      }
    } else {
      for(d=0; d<DIM; d++) {
	while(cg[d] < 0) {
	  cg[d] += box[d][d];
	  for(k=k0; (k<k1); k++) 
	    pos[cga[k]][d] += box[d][d];
	}
	while(cg[d] >= box[d][d]) {
	  cg[d] -= box[d][d];
	  for(k=k0; (k<k1); k++) 
	    pos[cga[k]][d] -= box[d][d];
	}
	cg_cm[icg][d] = cg[d];
      }
    }
#ifdef DEBUG_PBC
    for(d=0; (d<DIM); d++) {
      if ((cg_cm[icg][d] < 0) || (cg_cm[icg][d] >= box[d][d]))
	fatal_error(0,"cg_cm[%d] = %15f  %15f  %15f\n"
		    "box = %15f  %15f  %15f\n",
		    icg,cg_cm[icg][XX],cg_cm[icg][YY],cg_cm[icg][ZZ],
		    box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
    }
#endif
  }
}

void put_atoms_in_box(int natoms,matrix box,rvec x[])
{
  int i,m;

  for(i=0; (i<natoms); i++)
    for(m=0; m < DIM; m++) {
      while (x[i][m] < 0) 
	x[i][m] += box[m][m];
      while (x[i][m] >= box[m][m]) 
	x[i][m] -= box[m][m];
    }
}
