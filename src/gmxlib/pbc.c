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
 * Good ROcking Metal Altar for Chronical Sinners
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
static bool bTrunc;
/* static bool bTriclinic; */
static real side,side_2,side3_4;
/* static real diag; */
static rvec gl_fbox,gl_hbox,gl_mhbox;

void init_pbc(matrix box,bool bTruncOct)
{
  int i;
/*   int j; */

  bInit   = TRUE;
  bTrunc  = bTruncOct;
  side    = box[0][0];
  side_2  = 0.5*side;
  side3_4 = 0.75*side;
  /*
  diag    = side_2*sqrt(3.0);
  bTriclinic=FALSE;
  for(i=0; (i<DIM); i++)
    for(j=0; (j<DIM); j++)
      if ((i != j) && (box[i][j] != 0.0))
	bTriclinic=TRUE;
	*/
  for(i=0; (i<DIM); i++) {
    gl_fbox[i]  =  box[i][i];
    gl_hbox[i]  =  gl_fbox[i]*0.5;
    gl_mhbox[i] = -gl_hbox[i];
  }
}

static void do_trunc(rvec x)
{
  int i;

  if (fabs(x[XX]-side_2)+fabs(x[YY]-side_2)+fabs(x[ZZ]-side_2) > side3_4)
    for(i=0; (i<DIM); i++)
      x[i] -= sign(side_2,x[i]);
}

void pbc_dx(matrix box,rvec x1, rvec x2, rvec dx)
{
  int i;
  
  if (!bInit) { 
    fprintf(stderr,"init-pbc has not been called!\n");
    exit(1);
  }
  rvec_sub(x1,x2,dx);
  for(i=0; (i<DIM); i++) {
    if (dx[i] > gl_hbox[i])
      dx[i] -= gl_fbox[i];
    else if (dx[i] <= gl_mhbox[i])
      dx[i] += gl_fbox[i];
  }
  if (bTrunc)
    do_trunc(dx);
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

bool image_tri(ivec xi,ivec xj,imatrix box,real rlong2,int *shift,real *r2)
{
  fprintf(stderr,"image_tri still empty! (not implemented)\n");
  exit(1);

  return FALSE; /* To make the compiler happy */
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
	if (bTrunc && (k!=0) && (l!=0) && (m!=0)) {
	  /* Truncated edges... */
	  shift_vec[n][XX]=k*side_2;
	  shift_vec[n][YY]=l*side_2;
	  shift_vec[n][ZZ]=m*side_2;
	}
	else if (!bTrunc || (((k+l+m) % 2)!=0) || ((k==0) && (l==0) && (m==0)))
	  for(d = 0; d < DIM; d++)
	    shift_vec[n][d]=k*box[XX][d] + l*box[YY][d] + m*box[ZZ][d];
      }
}

void put_atoms_in_box(FILE *log,int cg0,int cg1,bool bTruncOct,
		      matrix box,rvec box_size,t_block *cgs,
		      rvec pos[],rvec shift_vec[],rvec cg_cm[])
{
  int  icg,ai,k,k0,k1;
  real dx,dy,dz,cgx,cgy,cgz,nrcg,inv_ncg;
  real bx,by,bz,hbx,hby,hbz,tx,ty,tz;
  real binv_x,binv_y,binv_z;
  atom_id *cga,*cgindex;
  
#define nint(x) ((x >= 0) ? ((int)(x+0.5)) : ((int)(x-0.5)))
  
#ifdef DEBUG
  fprintf(log,"Putting cgs %d to %d in box\n",cg0,cg1);
#endif
  bx      = box[XX][XX];
  by      = box[YY][YY];
  bz      = box[ZZ][ZZ];
  hbx     = 0.5*bx;
  hby     = 0.5*by;
  hbz     = 0.5*bz;
  binv_x  = 1.0/bx;
  binv_y  = 1.0/by;
  binv_z  = 1.0/bz;
  cga     = cgs->a;
  cgindex = cgs->index;
  
  for(icg=cg0; (icg<cg1); icg++) {
    /* First compute the center of geometry for this charge group */
    k0      = cgindex[icg];
    k1      = cgindex[icg+1];
    nrcg    = k1-k0;
    inv_ncg = 1.0/nrcg;
    
    cgx=cgy=cgz=0;
    for(k=k0; (k<k1); k++)  {
      ai=cga[k];
      
      cgx += inv_ncg*pos[ai][XX];
      cgy += inv_ncg*pos[ai][YY];
      cgz += inv_ncg*pos[ai][ZZ];
    }
    
    /* Now check pbc for this cg */
    tx = (cgx-hbx)*binv_x;
    ty = (cgy-hby)*binv_y;
    tz = (cgz-hbz)*binv_z;
    dx = bx*nint(tx);
    dy = by*nint(ty);
    dz = bz*nint(tz);

    /* We now have a shift vector for this CG 
     * let's add it to all the atoms and to the
     * cg_center of geometry
     */
    for(k=k0; (k<k1); k++) {
      ai=cga[k];
      
      pos[ai][XX] -= dx;
      pos[ai][YY] -= dy;
      pos[ai][ZZ] -= dz;
    }
    cg_cm[icg][XX] = cgx-dx;
    cg_cm[icg][YY] = cgy-dy;
    cg_cm[icg][ZZ] = cgz-dz;
  }
}

void put_all_atoms_in_box(int natoms,matrix box,rvec x[])
{
  int i,m;
 
  for(i=0; (i<natoms); i++)
    for(m=0; m < DIM; m++) {
      while (x[i][m] < 0) 
	x[i][m] += box[m][m];
      while (x[i][m] > box[m][m]) 
	x[i][m] -= box[m][m];
    }
}
