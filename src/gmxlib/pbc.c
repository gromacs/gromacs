/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
static ivec *tric_shift=NULL;
static matrix gl_box;
static int ntric_vec;
static real sure_dist2;

char *check_box(matrix box)
{
  char *ptr;

  if ((box[XX][YY]!=0) || (box[XX][ZZ]!=0) || (box[YY][ZZ]!=0))
    ptr = "Only triclinic boxes with the first vector parallel to the x-axis and the second vector in the xy-plane are supported.";
#ifdef ALLOW_OFFDIAG_LT_HALFDIAG 
  else if ((fabs(box[YY][XX])+fabs(box[ZZ][XX]) > 2*BOX_MARGIN*box[XX][XX]) ||
	   (fabs(box[ZZ][YY]) > 2*BOX_MARGIN*box[YY][YY]))
#else
  else if ((fabs(box[YY][XX]) > BOX_MARGIN*box[XX][XX]) ||
	   (fabs(box[ZZ][XX]) > BOX_MARGIN*box[XX][XX]) ||
	   (fabs(box[ZZ][YY]) > BOX_MARGIN*box[YY][YY]))
#endif
    ptr = "Triclinic box is too skewed.";
  else
    ptr = NULL;
  
  return ptr;
}

void init_pbc(matrix box)
{
  static int nalloc=0;
  int  i,j,k,d;
  real diagonal2;
  rvec try;
  char *ptr;

  ptr = check_box(box);
  if (ptr) {
    fprintf(stderr,   "Warning: %s\n",ptr);
    pr_rvecs(stderr,0,"         Box",box,DIM);
    fprintf(stderr,   "         Can not fix pbc.\n");
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
      sure_dist2 = 0.25*min(sqr(box[XX][XX]),
			    min(sqr(box[YY][YY]),sqr(box[ZZ][ZZ])));
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
		  srenew(tric_shift,nalloc);
		}
	      copy_rvec(try,tric_vec[ntric_vec]);
	      tric_shift[ntric_vec][XX] = i;
	      tric_shift[ntric_vec][YY] = j;
	      tric_shift[ntric_vec][ZZ] = k;
	      ntric_vec++;
	      }
	    }
    }
  }
  bInit   = TRUE;
}

int pbc_dx(const rvec x1, const rvec x2, rvec dx)
{
  int i,j;
  rvec dx_start,try;
  real d2min,d2try;
  ivec ishift;

  if (!bInit)
    fatal_error(0,"pbc_dx called before init_pbc");
  rvec_sub(x1,x2,dx);
  clear_ivec(ishift);
  if (bSupported) {
    if (bTriclinic) {
      for(i=DIM-1; i>=0; i--)
	if (dx[i] > gl_hbox[i]) {
	  for (j=i; j>=0; j--)
	    dx[j] -= gl_box[i][j];
	  ishift[i]--;
	} else if (dx[i] <= gl_mhbox[i]) {
	  for (j=i; j>=0; j--)
	    dx[j] += gl_box[i][j];
	  ishift[i]++;
	}
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
	  ivec_inc(ishift,tric_shift[i]);
	  d2min = d2try;
	}
	i++;
      }
    } else {
      for(i=0; i<DIM; i++)
	if (dx[i] > gl_hbox[i]) {
	  dx[i] -= gl_fbox[i];
	  ishift[i]--;
	} else if (dx[i] <= gl_mhbox[i]) {
	  dx[i] += gl_fbox[i];
	  ishift[i]++;
	}
    }
  }

  return IVEC2IS(ishift);
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

void calc_shifts(matrix box,rvec box_size,rvec shift_vec[])
{
  int k,l,m,d,n,test;
  
  init_pbc(box);
  for (m=0; (m<DIM); m++)
    box_size[m]=box[m][m];
  
  n=0;
  
  for(m = -D_BOX_Z; m <= D_BOX_Z; m++)
    for(l = -D_BOX_Y; l <= D_BOX_Y; l++) 
      for(k = -D_BOX_X; k <= D_BOX_X; k++,n++) {
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

void put_charge_groups_in_box(FILE *log,int cg0,int cg1,
			      matrix box,rvec box_size,t_block *cgs,
			      rvec pos[],rvec cg_cm[])
			      
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

void calc_box_center(matrix box,rvec box_center)
{
  int d,m;

  clear_rvec(box_center);
  for(m=0; (m<DIM); m++)  
    for(d=0; d<DIM; d++)
      box_center[d] += 0.5*box[m][d];
}

void calc_triclinic_images(matrix box,rvec img[])
{
  int i;

  /* Calculate 3 adjacent images in the xy-plane */
  copy_rvec(box[0],img[0]);
  copy_rvec(box[1],img[1]);
  if (img[1][XX] < 0)
    svmul(-1,img[1],img[1]);
  rvec_sub(img[1],img[0],img[2]);

  /* Get the next 3 in the xy-plane as mirror images */
  for(i=0; i<3; i++)
    svmul(-1,img[i],img[3+i]);

  /* Calculate the first 4 out of xy-plane images */
  copy_rvec(box[2],img[6]);
  if (img[6][XX] < 0)
    svmul(-1,img[6],img[6]);
  for(i=0; i<3; i++)
    rvec_add(img[6],img[i+1],img[7+i]);

  /* Mirror the last 4 from the previous in opposite rotation */
  for(i=0; i<4; i++)
    svmul(-1,img[6 + (2+i) % 4],img[10+i]);
}

void calc_compact_unitcell_vertices(matrix box,rvec vert[])
{
  rvec img[NTRICIMG],box_center;
  int n,i,j,tmp[4],d;

  calc_triclinic_images(box,img);

  n=0;
  for(i=2; i<=5; i+=3) {
    tmp[0] = i-1;
    if (i==2)
      tmp[1] = 8;
    else 
      tmp[1] = 6;
    tmp[2] = (i+1) % 6;
    tmp[3] = tmp[1]+4;
    for(j=0; j<4; j++) {
      for(d=0; d<DIM; d++)
	vert[n][d] = img[i][d]+img[tmp[j]][d]+img[tmp[(j+1)%4]][d];
      n++;
    }
  }
  for(i=7; i<=13; i+=6) {
    tmp[0] = (i-7)/2;
    tmp[1] = tmp[0]+1;
    if (i==7)
      tmp[2] = 8;
    else
      tmp[2] = 10;
    tmp[3] = i-1;
    for(j=0; j<4; j++) {
      for(d=0; d<DIM; d++)
	vert[n][d] = img[i][d]+img[tmp[j]][d]+img[tmp[(j+1)%4]][d];
      n++;
    }
  }
  for(i=9; i<=11; i+=2) {
    if (i==9)
      tmp[0] = 3;
    else
      tmp[0] = 0;
    tmp[1] = tmp[0]+1;
    if (i==9)
      tmp[2] = 6;
    else
      tmp[2] = 12;
    tmp[3] = i-1;
    for(j=0; j<4; j++) {
      for(d=0; d<DIM; d++)
	vert[n][d] = img[i][d]+img[tmp[j]][d]+img[tmp[(j+1)%4]][d];
      n++;
    }
  }

  calc_box_center(box,box_center);
  for(i=0; i<NCUCVERT; i++)
    for(d=0; d<DIM; d++)
      vert[i][d] = vert[i][d]*0.25+box_center[d];
}

int *compact_unitcell_edges()
{
  /* this is an index in vert[] (see calc_box_vertices) */
  static int edge[NCUCEDGE*2];
  static int hexcon[24] = { 0,9, 1,19, 2,15, 3,21, 
			    4,17, 5,11, 6,23, 7,13,
			    8,20, 10,18, 12,16, 14,22 };
  int e,i,j;
  bool bFirst = TRUE;

  if (bFirst) {
    e = 0;
    for(i=0; i<6; i++)
      for(j=0; j<4; j++) {
	edge[e++] = 4*i + j;
	edge[e++] = 4*i + (j+1) % 4;
      }
    for(i=0; i<12*2; i++)
      edge[e++] = hexcon[i];
    
    bFirst = FALSE;
  }

  return edge;
}

void put_atoms_in_box(matrix box,int natoms,rvec x[])
{
  int i,m,d;

  for(i=0; (i<natoms); i++)
    for(m=DIM-1; m>=0; m--) {
      while (x[i][m] < 0) 
	for(d=0; d<=m; d++)
	  x[i][d] += box[m][d];
      while (x[i][m] >= box[m][m])
	for(d=0; d<=m; d++)
	  x[i][d] -= box[m][d];
    }
}

void put_atoms_in_triclinic_unitcell(matrix box,int natoms,rvec x[])
{
  rvec   box_center,shift_center;
  real   shm01,shm02,shm12,shift;
  int    i,m,d;
  
  calc_box_center(box,box_center);
  
  /* The product of matrix shm with a coordinate gives the shift vector
     which is required determine the periodic cell position */
  shm01 = box[1][0]/box[1][1];
  shm02 = (box[1][1]*box[2][0] - box[2][1]*box[1][0])/(box[1][1]*box[2][2]);
  shm12 = box[2][1]/box[2][2];

  clear_rvec(shift_center);
  for(d=0; d<DIM; d++)
    rvec_inc(shift_center,box[d]);
  svmul(0.5,shift_center,shift_center);
  rvec_sub(box_center,shift_center,shift_center);

  shift_center[0] = shm01*shift_center[1] + shm02*shift_center[2];
  shift_center[1] = shm12*shift_center[2];
  shift_center[2] = 0;

  for(i=0; (i<natoms); i++)
    for(m=DIM-1; m>=0; m--) {
      shift = shift_center[m];
      if (m == 0) {
	shift += shm01*x[i][1] + shm02*x[i][2];
      } else if (m == 1) {
	shift += shm12*x[i][2];
      }
      while (x[i][m]-shift < 0)
	for(d=0; d<=m; d++)
	  x[i][d] += box[m][d];
      while (x[i][m]-shift >= box[m][m])
	for(d=0; d<=m; d++)
	  x[i][d] -= box[m][d];
    }
}

void put_atoms_in_compact_unitcell(matrix box,int natoms,rvec x[])
{
  rvec box_center,dx;
  int  i;

  init_pbc(box);
  calc_box_center(box,box_center);
  for(i=0; i<natoms; i++) {
    pbc_dx(x[i],box_center,dx);
    rvec_add(box_center,dx,x[i]);
  }
}

