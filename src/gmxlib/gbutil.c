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
static char *SRCID_gbutil_c = "$Id$";

#include <math.h>
#include "macros.h"
#include "vec.h"
#include "fatal.h"
#include "gstat.h"
#include "pbc.h"

real dist2(rvec x,rvec y,matrix box)
{
  rvec dx;
  
  pbc_dx(box,x,y,dx);
  
  return iprod(dx,dx);
}

real distance_to_z(rvec x)
{
  return (sqr(x[XX])+sqr(x[YY]));
} /*distance_to_z()*/

int in_box(int NTB,matrix box,rvec x)
{
  switch (NTB) {
  /* cubic box*/
  case 0: 
    if (
	( (x[0]>=0.0) && (x[0]<=box[0][0]) ) &&
	( (x[1]>=0.0) && (x[1]<=box[1][1]) ) &&
	( (x[2]>=0.0) && (x[2]<=box[2][2]) )
	)
      return 0;
    else
      return 1;
      
  /*truncated octahedron*/
  case 1: 
    if (
	((fabs(box[XX][XX]/2-x[XX]) <= box[XX][XX]/2) &&
	 (fabs(box[YY][YY]/2-x[YY]) <= box[YY][YY]/2) &&
	 (fabs(box[ZZ][ZZ]/2-x[ZZ]) <= box[ZZ][ZZ]/2))
	&&
	(fabs(x[XX]-box[XX][XX]/2)+fabs(x[YY]-box[YY][YY]/2)+fabs(x[ZZ]-box[ZZ][ZZ]/2) 
	 <3*box[XX][XX]/4)
	)
      return 0;
    else 
      return 1;
  default: 
    fatal_error(0,"Illegal boxtype\nProgram terminated");
  }
  return 1;
}/*in_box()*/


void rotate_conf(int natom,rvec *x,rvec *v,real alfa, real beta,real gamma)
{
  int  i,m;
  rvec x_old,v_old;
  
  for (i=0; (i<natom); i++) { 
    /*copy x[i] to x_old*/
    for(m=0;(m<DIM);m++) {
      x_old[m]=x[i][m];
      v_old[m]=v[i][m];
    }
    
    /*calculate new x[i] by rotation alfa around the x-axis*/
    x[i][XX]=x_old[XX];
    v[i][XX]=v_old[XX];
    x[i][YY]=x_old[YY]*cos(alfa)+x_old[ZZ]*sin(alfa);
    v[i][YY]=v_old[YY]*cos(alfa)+v_old[ZZ]*sin(alfa);
    x[i][ZZ]=x_old[ZZ]*cos(alfa)-x_old[YY]*sin(alfa);
    v[i][ZZ]=v_old[ZZ]*cos(alfa)-v_old[YY]*sin(alfa);
    
    /*copy x[i] to x_old*/ 
    for(m=0;(m<DIM);m++) {
      x_old[m]=x[i][m];
      v_old[m]=v[i][m];
    }
    
    /*calculate new x[i] by rotation beta around the y-axis*/
    x[i][XX]=x_old[XX]*cos(beta)-x_old[ZZ]*sin(beta);
    x[i][YY]=x_old[YY];
    x[i][ZZ]=x_old[XX]*sin(beta)+x_old[ZZ]*cos(beta);
    v[i][XX]=v_old[XX]*cos(beta)-v_old[ZZ]*sin(beta);
    v[i][YY]=v_old[YY];
    v[i][ZZ]=v_old[XX]*sin(beta)+v_old[ZZ]*cos(beta);
    
    /*copy x[i] to x_old*/
    for(m=0;(m<DIM);m++) {
      x_old[m]=x[i][m];
      v_old[m]=v[i][m];
    }
    
    /*calculate new x[i] by rotation gamma around the z-axis*/
    x[i][XX]=x_old[XX]*cos(gamma)+x_old[YY]*sin(gamma);
    x[i][YY]=x_old[YY]*cos(gamma)-x_old[XX]*sin(gamma);
    x[i][ZZ]=x_old[ZZ];
    v[i][XX]=v_old[XX]*cos(gamma)+v_old[YY]*sin(gamma);
    v[i][YY]=v_old[YY]*cos(gamma)-v_old[XX]*sin(gamma);
    v[i][ZZ]=v_old[ZZ];
  }
}/*rotate_conf()*/


void orient(int natom,rvec *x,rvec *v, rvec angle,matrix box)
{
  real longest,rij,rzi;
  int  i,j,m,max_i=0,max_j=0;
  rvec origin;
  int  temp;
  real alfa=0,beta=0,gamma=0;
  
  /*first i am going to look for the longest atom-atom distance*/
  longest=dist2(x[0],x[1],box);
  i=0;
  j=1;
  for (i=0;(i<natom);i++) {
    for (j=0;(j<natom);j++) {
      rij=dist2(x[i],x[j],box);
      if (rij>longest) {
	max_i=i;
	max_j=j;
	longest=rij;
      }
    }
  }
  /* first check if x[max_i]<x[max_j] else swap*/
  if (x[max_i][2]>x[max_j][2]) {
    temp=max_i;
    max_i=max_j;
    max_j=temp;
  }
  
  /*set the origin to x[i]*/
  for(m=0;(m<DIM);m++) 
    origin[m]=x[max_i][m];
  for(i=0;(i<natom);i++)
    for(m=0;(m<DIM);m++)
      x[i][m]-=origin[m];
      
  /* calculate the rotation angles alfa(x_axis) and beta(y_axis)
   * the rotation angles must be calculated clockwise looking 
   * along the rotation axis to the origin*
   * alfa (x-axis)
   */
  alfa=atan(x[max_j][ZZ]/x[max_j][YY])-M_PI_2;
  beta=M_PI_2-atan(x[max_j][ZZ]/x[max_j][XX]);
  rotate_conf(natom,x,v,alfa,beta,gamma);
  
  /* now search the longest distance for rotation along the z_axis */
  longest=distance_to_z(x[0]);
  max_i=0;
  for (i=1;(i<natom);i++) {
    rzi=distance_to_z(x[i]);
    if (rzi>longest) {
      longest = rzi;
      max_i=i;
    }
  }
  gamma=atan(x[max_i][YY]/x[max_i][XX])-M_PI_2;
  rotate_conf(natom,x,v,0,0,gamma);
  angle[0]=alfa;
  angle[1]=beta;
  angle[2]=gamma;
} /*orient()*/


void genconf(t_atoms *atoms,rvec *x,rvec *v,real *r,matrix box,ivec n_box)
{
  int     i,ix,iy,iz,m,j,imol=0,offset;
  rvec    delta;
  int     nmol=1;
  
  for (i=0;i<DIM;i++)
    nmol*=n_box[i];
    
  /*print message*/
  fprintf(stderr,"Generating configuration\n");
  for(ix=0; (ix < n_box[0]); ix++) {
    delta[XX]=ix*box[XX][XX];
    for(iy=0; (iy < n_box[1]); iy++) {
      delta[YY]=iy*box[YY][YY];
      for(iz=0; (iz < n_box[2]); iz++) {
	delta[ZZ]=iz*box[ZZ][ZZ];
	offset=imol*atoms->nr;
	for (i=0;(i < atoms->nr);i++) {
	  for (m=0;(m < DIM);m++) {
	    x[offset+i][m]=delta[m]+x[i][m];
	    v[offset+i][m]=v[i][m];
	  }
	  r[offset+i]=r[i];
        }
	imol++;
      }
    }
  }
  for (i=1;(i<nmol);i++) {
    int offs    = i*atoms->nr;
    int offsres = i*atoms->nres;
    for (j=0;(j<atoms->nr);j++) {
      atoms->atomname[offs+j]  = atoms->atomname[j];
      atoms->atom[offs+j].resnr = atoms->atom[j].resnr+offsres;
      atoms->resname[atoms->atom[offs+j].resnr]=
	atoms->resname[atoms->atom[j].resnr];
    }
  }
  atoms->nr*=nmol;
  atoms->nres*=nmol;
  for(i=0;(i<DIM);i++) {
    box[XX][i]*=n_box[0];
    box[YY][i]*=n_box[1];
    box[ZZ][i]*=n_box[2];
  }
} /*genconf()*/

/*gen_box() generates a box around a configuration*/
void gen_box(int NTB,int natoms,rvec *x, matrix box,rvec box_space,
	     bool bCenter)
{
  int i,m;
  rvec xmin, xmax;
  real max_box;
  
  /*calculate minimum and maximum x[0..DIM-1]*/
  for (m=0;(m<DIM);m++)
    xmin[m]=xmax[m]=x[0][m];
  for (i=1;(i < natoms); i++) 
    for (m=0;m<DIM;m++) {
      xmin[m]=min(xmin[m],x[i][m]);
      xmax[m]=max(xmax[m],x[i][m]);
    }
    
  /*calculate the new box sizes for cubic and octahedral ...*/
  for (m=0; (m<DIM);m++)
    box[m][m]=xmax[m]-xmin[m]+2*box_space[m];
 
  /*calculate the box size if NTB=1 (truncated octahedron)*/
  if (NTB==1) {
    max_box=box[0][0];
    for(m=0;(m<DIM);m++)
      max_box=max(max_box,box[m][m]); 
    for (m=0;(m<DIM);m++)
      box[m][m]=max_box;
  }
  
  /*move the molecule to the center of the box*/
  if (bCenter)
    for(i=0;(i<natoms);i++)
      for (m=0;(m<DIM);m++) {
	x[i][m]+=0.5*(box[m][m]-xmin[m]-xmax[m]);
      }


#ifdef DEBUG 
  /* print data to check this */
  print_stat(x,natoms,box);
#endif
}/*gen_box()*/

