/*
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "statutil.h"
#include "macros.h"
#include "tpxio.h"
#include "smalloc.h"
#include "physics.h"
#include "vec.h"
#include "gstat.h"
#include "nrjac.h"
#include "copyrite.h"
#include "index.h"
#include "gmx_ana.h"


#define NM2ANG 10
#define TOLERANCE 1.0E-8

#define e2d(x) ENM2DEBYE*(x)
#define delta(a,b) (( a == b ) ? 1.0 : 0.0)

#define NDIM 3          /* We will be using a numerical recipes routine */

static char dim[DIM+1] = "XYZ";

typedef real        	tensor3[DIM][DIM][DIM];      /* 3 rank tensor */
typedef real        	tensor4[DIM][DIM][DIM][DIM]; /* 4 rank tensor */


void pr_coord(int k0,int k1,atom_id index[],rvec x[],char *msg) 
{
  int k,kk;

  fprintf(stdout,"Coordinates in nm (%s)\n",msg);
  for(k=k0; (k<k1); k++) {
    kk=index[k];
    fprintf(stdout,"Atom %d, %15.10f %15.10f %15.10f\n",
	    kk,x[kk][XX],x[kk][YY],x[kk][ZZ]);
  }
  fprintf(stdout,"\n");
}

static void clear_tensor3(tensor3 a)
{
  int i,j,k;
  const real nul=0.0;
  
  for(i=0; (i<DIM); i++)
    for(j=0; (j<DIM); j++)
      for(k=0; (k<DIM); k++)
	a[i][j][k]=nul;
}

static void clear_tensor4(tensor4 a)
{
  int i,j,k,l;
  const real nul=0.0;
  
  for(i=0; (i<DIM); i++)
    for(j=0; (j<DIM); j++)
      for(k=0; (k<DIM); k++)
	for(l=0; (l<DIM); l++)
	  a[i][j][k][l]=nul;
}

void rotate_mol(int k0,int k1,atom_id index[],rvec x[],matrix trans)
{
  real   xt,yt,zt;
  int    k,kk;

  for(k=k0; (k<k1); k++) {
    kk=index[k];
    xt=x[kk][XX];
    yt=x[kk][YY];
    zt=x[kk][ZZ];
    x[kk][XX]=trans[XX][XX]*xt+trans[XX][YY]*yt+trans[XX][ZZ]*zt;
    x[kk][YY]=trans[YY][XX]*xt+trans[YY][YY]*yt+trans[YY][ZZ]*zt;
    x[kk][ZZ]=trans[ZZ][XX]*xt+trans[ZZ][YY]*yt+trans[ZZ][ZZ]*zt;
  }
}

/* the following routines are heavily inspired by the Gaussian 94 source
 * code 
 */

/*
     Make the rotation matrix for angle Theta counterclockwise about          
     axis IXYZ.                                                              
*/

void make_rot_mat(int axis,real theta,matrix t_mat){
  
  ivec i;
  real s,c;

  
  i[XX]=axis + 1;
  i[YY]=1 + i[XX] % 3;
  i[ZZ]=1 + i[YY] % 3;

  i[XX]-=1;
  i[YY]-=1;
  i[ZZ]-=1;

  s=sin(theta);
  c=cos(theta);
  t_mat[i[XX]][i[XX]]=1.0;
  t_mat[i[XX]][i[YY]]=0.0;
  t_mat[i[XX]][i[ZZ]]=0.0;
  t_mat[i[YY]][i[XX]]=0.0;
  t_mat[i[YY]][i[YY]]=c;
  t_mat[i[YY]][i[ZZ]]=s;
  t_mat[i[ZZ]][i[XX]]=0.0;
  t_mat[i[ZZ]][i[YY]]=-s;
  t_mat[i[ZZ]][i[ZZ]]=c;
}

gmx_bool test_linear_mol(rvec d)
{
  /* d is sorted in descending order */
  if ( (d[ZZ] < TOLERANCE) && (d[XX]-d[YY]) < TOLERANCE ) {
    return TRUE;
  } else 
    return FALSE;
}

/* Returns the third moment of charge along an axis */
real test_qmom3(int k0,int k1,atom_id index[],t_atom atom[],rvec x[],int axis){

  int k,kk;
  real xcq,q;

  xcq=0.0;
  for(k=k0; (k<k1); k++) {
    kk=index[k];
    q=fabs(atom[kk].q);
    xcq+=q*x[kk][axis]*x[kk][axis]*x[kk][axis];
  }
  
  return xcq;
}

/* Returns the second moment of mass along an axis */
real test_mmom2(int k0,int k1,atom_id index[],t_atom atom[],rvec x[],int axis){

  int k,kk;
  real xcm,m;

  xcm=0.0;
  for(k=k0; (k<k1); k++) {
    kk=index[k];
    m=atom[kk].m;
    xcm+=m*x[kk][axis]*x[kk][axis];
  }
  
  return xcm;
}

real calc_xcm_mol(int k0,int k1,atom_id index[],t_atom atom[],rvec x[],
		  rvec xcm)
{
  int  k,kk,m;
  real m0,tm;

  /* Compute the center of mass */
  clear_rvec(xcm);
  tm=0.0;
  for(k=k0; (k<k1); k++) {
    kk=index[k];
    m0=atom[kk].m;
    tm+=m0;
    for(m=0; (m<DIM); m++)
      xcm[m]+=m0*x[kk][m];
  }
  for(m=0; (m<DIM); m++)
    xcm[m]/=tm;

  /* And make it the origin */
  for(k=k0; (k<k1); k++) {
    kk=index[k];
    rvec_dec(x[kk],xcm);
  }
  
  return tm;
}

/* Retruns the center of charge */
real calc_xcq_mol(int k0,int k1,atom_id index[],t_atom atom[],rvec x[],
		  rvec xcq)
{
  int  k,kk,m;
  real q0,tq;

  clear_rvec(xcq);
  tq=0.0;
  for(k=k0; (k<k1); k++) {
    kk=index[k];
    q0=fabs(atom[kk].q);
    tq+=q0;
    fprintf(stdout,"tq: %f, q0: %f\n",tq,q0);
    for(m=0; (m<DIM); m++)
      xcq[m]+=q0*x[kk][m];
  }

  for(m=0; (m<DIM); m++)
    xcq[m]/=tq;
  /*
  for(k=k0; (k<k1); k++) {
    kk=index[k];
    rvec_dec(x[kk],xcq);
  }
  */
  return tq;
}

/* Returns in m1 the dipole moment */
void mol_M1(int n0,int n1,atom_id ma[],rvec x[],t_atom atom[],rvec m1)
{
  int  m,n,nn;
  real q;

  clear_rvec(m1);
  for(n=n0; (n<n1); n++) {
    nn = ma[n];
    q  = e2d(atom[nn].q);
    for(m=0; (m<DIM); m++)
      m1[m] += q*x[nn][m];
  }
}

/* returns in m2 the quadrupole moment */
void mol_M2(int n0,int n1,atom_id ma[],rvec x[],t_atom atom[],tensor m2)
{
  int  n,nn,i,j;
  real q,r2;

  clear_mat(m2);
  for(n=n0; (n<n1); n++) {
    nn = ma[n];
    q  = e2d(atom[nn].q);
    r2 = norm2(x[nn]);
    for(i=0; (i<DIM); i++)
      for(j=0; (j<DIM); j++)
	m2[i][j] += 0.5*q*(3.0*x[nn][i]*x[nn][j] - r2*delta(i,j))*NM2ANG;
  }
}

/* Returns in m3 the octopole moment */
void mol_M3(int n0,int n1,atom_id ma[],rvec x[],t_atom atom[],tensor3 m3)
{
  int  i,j,k,n,nn;
  real q,r2;

  clear_tensor3(m3);
  for(n=n0; (n<n1); n++) {
    nn = ma[n];
    q  = e2d(atom[nn].q);
    r2 = norm2(x[nn]);
    for(i=0; (i<DIM); i++)
      for(j=0; (j<DIM); j++)
	for(k=0; (k<DIM); k++)
	  m3[i][j][k] += 
	    0.5*q*(5.0*x[nn][i]*x[nn][j]*x[nn][k] 
		   - r2*(x[nn][i]*delta(j,k) + 
			 x[nn][j]*delta(k,i) +
			 x[nn][k]*delta(i,j)))*NM2ANG*NM2ANG;
  }
}

/* Returns in m4 the hexadecapole moment */
void mol_M4(int n0,int n1,atom_id ma[],rvec x[],t_atom atom[],tensor4 m4)
{
  int  i,j,k,l,n,nn;
  real q,r2;

  clear_tensor4(m4);
  for(n=n0; (n<n1); n++) {
    nn = ma[n];
    q  = e2d(atom[nn].q);
    r2 = norm2(x[nn]);
    for(i=0; (i<DIM); i++)
      for(j=0; (j<DIM); j++)
	for(k=0; (k<DIM); k++)
	  for(l=0; (l<DIM); l++)
	    m4[i][j][k][l] += 
	      0.125*q*(35.0*x[nn][i]*x[nn][j]*x[nn][k]*x[nn][l] 
		     - 5.0*r2*(x[nn][i]*x[nn][j]*delta(k,l) + 
			       x[nn][i]*x[nn][k]*delta(j,l) +
			       x[nn][i]*x[nn][l]*delta(j,k) +
			       x[nn][j]*x[nn][k]*delta(i,l) +
			       x[nn][j]*x[nn][l]*delta(i,k) +
			       x[nn][k]*x[nn][l]*delta(i,j)) 
		      + r2*r2*(delta(i,j)*delta(k,l) + 
			       delta(i,k)*delta(j,l) +
			       delta(i,l)*delta(j,k)))*NM2ANG*NM2ANG*NM2ANG;
  }
}

/* Print the dipole moment components and the total dipole moment */
void pr_M1(FILE *fp,char *msg,int mol,rvec m1,real time)
{
  int  i;
  real m1_tot;

  fprintf(fp,"Molecule: %d @ t= %f ps\n",mol,time);
  
  m1_tot = sqrt(m1[XX]*m1[XX]+m1[YY]*m1[YY]+m1[ZZ]*m1[ZZ]);
  
  fprintf(stdout,"Dipole Moment %s(Debye):\n",msg);
  fprintf(stdout,"X= %10.5f Y= %10.5f Z= %10.5f Tot= %10.5f\n",
	  m1[XX],m1[YY],m1[ZZ],m1_tot);
}

/* Print the quadrupole moment components */
void pr_M2(FILE *fp,char *msg,tensor m2,gmx_bool bFull)
{
  int i,j;

  fprintf(fp,"Quadrupole Moment %s(Debye-Ang):\n",msg);
  if (!bFull) {
    fprintf(fp,"XX= %10.5f YY= %10.5f ZZ= %10.5f\n",
	    m2[XX][XX],m2[YY][YY],m2[ZZ][ZZ]);
    fprintf(fp,"XY= %10.5f XZ= %10.5f YZ= %10.5f\n",
	    m2[XX][YY],m2[XX][ZZ],m2[YY][ZZ]);
  }
  else {
    for(i=0; (i<DIM); i++) {
      for(j=0; (j<DIM); j++)
	fprintf(fp,"  %c%c= %10.4f",dim[i],dim[j],m2[i][j]);
      fprintf(fp,"\n");
    }
  }
}

/* Print the octopole moment components */
void pr_M3(FILE *fp,char *msg,tensor3 m3,gmx_bool bFull)
{
  int i,j,k;

  fprintf(fp,"Octopole Moment %s(Debye-Ang^2):\n",msg);
  if (!bFull) {
    fprintf(fp,"XXX= %10.5f YYY= %10.5f ZZZ= %10.5f XYY= %10.5f\n",
	    m3[XX][XX][XX],m3[YY][YY][YY],m3[ZZ][ZZ][ZZ],m3[XX][YY][YY]);
    fprintf(fp,"XXY= %10.5f XXZ= %10.5f XZZ= %10.5f YZZ= %10.5f\n",
	    m3[XX][XX][YY],m3[XX][XX][ZZ],m3[XX][ZZ][ZZ],m3[YY][ZZ][ZZ]);
    fprintf(fp,"YYZ= %10.5f XYZ= %10.5f\n",
	    m3[YY][YY][ZZ],m3[XX][YY][ZZ]);
  }
  else {
    for(i=0; (i<DIM); i++) {
      for(j=0; (j<DIM); j++) {
	for(k=0; (k<DIM); k++)
	  fprintf(fp,"  %c%c%c= %10.4f",dim[i],dim[j],dim[k],m3[i][j][k]);
	fprintf(fp,"\n");
      }
    }
  }
}

/* Print the hexadecapole moment components */
void pr_M4(FILE *fp,char *msg,tensor4 m4,gmx_bool bFull)
{
  int i,j,k,l;

  fprintf(fp,"Hexadecapole Moment %s(Debye-Ang^3):\n",msg);
  if (!bFull) {
    fprintf(fp,"XXXX= %10.5f YYYY= %10.5f ZZZZ= %10.5f XXXY= %10.5f\n",
	    m4[XX][XX][XX][XX],m4[YY][YY][YY][YY],
	    m4[ZZ][ZZ][ZZ][ZZ],m4[XX][XX][XX][YY]);
    fprintf(fp,"XXXZ= %10.5f YYYX= %10.5f YYYZ= %10.5f ZZZX= %10.5f\n",
	    m4[XX][XX][XX][ZZ],m4[YY][YY][YY][XX],
	    m4[YY][YY][YY][ZZ],m4[ZZ][ZZ][ZZ][XX]);
    fprintf(fp,"ZZZY= %10.5f XXYY= %10.5f XXZZ= %10.5f YYZZ= %10.5f\n",
	    m4[ZZ][ZZ][ZZ][YY],m4[XX][XX][YY][YY],
	    m4[XX][XX][ZZ][ZZ],m4[YY][YY][ZZ][ZZ]);
    fprintf(fp,"XXYZ= %10.5f YYXZ= %10.5f ZZXY= %10.5f\n\n",
	    m4[XX][XX][YY][ZZ],m4[YY][YY][XX][ZZ],m4[ZZ][ZZ][XX][YY]);
  }
  else {
    for(i=0; (i<DIM); i++) {
      for(j=0; (j<DIM); j++) {
	for(k=0; (k<DIM); k++) {
	  for(l=0; (l<DIM); l++)
	    fprintf(fp,"  %c%c%c%c = %10.4f",dim[i],dim[j],dim[k],dim[l],
		    m4[i][j][k][l]);
	  fprintf(fp,"\n");
	}
      }
    }
  }
}

/* Compute the inertia tensor and returns in trans a matrix which rotates
 * the molecules along the principal axes system */
void principal_comp_mol(int k0,int k1,atom_id index[],t_atom atom[],rvec x[],
			matrix trans,rvec d)
{
  int  i,j,ai,m,nrot;
  real mm,rx,ry,rz;
  double **inten,dd[NDIM],tvec[NDIM],**ev;
  real temp;
  
  snew(inten,NDIM);
  snew(ev,NDIM);
  for(i=0; (i<NDIM); i++) {
    snew(inten[i],NDIM);
    snew(ev[i],NDIM);
    dd[i]=0.0;
  }
    
  for(i=0; (i<NDIM); i++)
    for(m=0; (m<NDIM); m++)
      inten[i][m]=0;
 
  for(i=k0; (i<k1); i++) {
    ai=index[i];
    mm=atom[ai].m;
    rx=x[ai][XX];
    ry=x[ai][YY];
    rz=x[ai][ZZ];
    inten[0][0]+=mm*(sqr(ry)+sqr(rz));
    inten[1][1]+=mm*(sqr(rx)+sqr(rz));
    inten[2][2]+=mm*(sqr(rx)+sqr(ry));
    inten[1][0]-=mm*(ry*rx);
    inten[2][0]-=mm*(rx*rz);
    inten[2][1]-=mm*(rz*ry);
  }
  inten[0][1]=inten[1][0];
  inten[0][2]=inten[2][0];
  inten[1][2]=inten[2][1];
  
  /* Call numerical recipe routines */
  jacobi(inten,3,dd,ev,&nrot);
  
  /* Sort eigenvalues in descending order */
#define SWAPPER(i) 			\
  if (fabs(dd[i+1]) > fabs(dd[i])) {	\
    temp=dd[i];			\
    for(j=0; (j<NDIM); j++) tvec[j]=ev[j][i];\
    dd[i]=dd[i+1];			\
    for(j=0; (j<NDIM); j++) ev[j][i]=ev[j][i+1];		\
    dd[i+1]=temp;			\
    for(j=0; (j<NDIM); j++) ev[j][i+1]=tvec[j];			\
  }
  SWAPPER(0)
  SWAPPER(1)
  SWAPPER(0)
      
  for(i=0; (i<DIM); i++) {
    d[i]=dd[i];
    for(m=0; (m<DIM); m++)
      trans[i][m]=ev[m][i];
  }
    
  for(i=0; (i<NDIM); i++) {
    sfree(inten[i]);
    sfree(ev[i]);
  }
  sfree(inten);
  sfree(ev);
}


/* WARNING WARNING WARNING
 * This routine rotates a molecule (I have checked this for water, PvM)
 * in the standard orientation used for water by researchers in the field. 
 * This is different from the orientation used by Gray and Gubbins, 
 * so be careful, with molecules other than water */
void rot_mol_to_std_orient(int k0,int k1,atom_id index[],t_atom atom[],
			   rvec x[],matrix trans)
{
  int  i;
  rvec xcm,xcq,d;
  matrix r_mat;

  clear_rvec(xcm);

  /* Compute the center of mass of the molecule and make it the origin */
  calc_xcm_mol(k0,k1,index,atom,x,xcm);
  
  /* Compute the inertia moment tensor of a molecule */
  principal_comp_mol(k0,k1,index,atom,x,trans,d);

  /* Rotate molecule to align with principal axes */
  rotate_mol(k0,k1,index,x,trans); 

  
  /* If one of the moments is zero and the other two are equal, the 
   * molecule is linear
   */
  
  if (test_linear_mol(d)) {
    fprintf(stdout,"This molecule is linear\n");
  } else {
    fprintf(stdout,"This molecule is not linear\n");

make_rot_mat(ZZ,-0.5*M_PI,r_mat);
rotate_mol(k0,k1,index,x,r_mat);



    /* Now check if the center of charge now lies on the Z-axis 
     * If not, rotate molecule so that it does.
     */
    for(i=0; (i<DIM); i++) {
      xcq[i]=test_qmom3(k0,k1,index,atom,x,i);
    }

    if ((fabs(xcq[ZZ]) - TOLERANCE) < 0.0) {
      xcq[ZZ]=0.0;
    } else {
#ifdef DEBUG
      fprintf(stdout,"Center of charge on Z-axis: %f\n",xcq[ZZ]);
#endif
      if (xcq[ZZ] > 0.0) {
	make_rot_mat(XX,M_PI,r_mat);
	rotate_mol(k0,k1,index,x,r_mat); 
      }
    }
    
    if ((fabs(xcq[XX]) - TOLERANCE) < 0.0) {
      xcq[XX]=0.0;
    } else {
#ifdef DEBUG
      fprintf(stdout,"Center of charge on X-axis: %f\n",xcq[XX]);
#endif
      if (xcq[XX] < 0.0) {
	make_rot_mat(YY,0.5*M_PI,r_mat);
	rotate_mol(k0,k1,index,x,r_mat); 
      } else {
	make_rot_mat(YY,-0.5*M_PI,r_mat);
	rotate_mol(k0,k1,index,x,r_mat); 
      }
    }
      
    if ((fabs(xcq[YY]) - TOLERANCE) < 0.0) {
      xcq[YY]=0.0;
    } else {
#ifdef DEBUG
      fprintf(stdout,"Center of charge on Y-axis: %f\n",xcq[YY]);
#endif
      if (xcq[YY] < 0.0) {
	make_rot_mat(XX,-0.5*M_PI,r_mat);
	rotate_mol(k0,k1,index,x,r_mat);
      } else {
	make_rot_mat(XX,0.5*M_PI,r_mat);
	rotate_mol(k0,k1,index,x,r_mat);
      }
    }

    /* Now check the trace of the inertia tensor.
     * I want the water molecule in the YZ-plane */
    for(i=0; (i<DIM); i++) {
      xcm[i]=test_mmom2(k0,k1,index,atom,x,i);
    }
#ifdef DEBUG
    fprintf(stdout,"xcm: %f %f %f\n",xcm[XX],xcm[YY],xcm[ZZ]);
#endif
    
    /* Check if X-component of inertia tensor is zero, if not 
     * rotate molecule
     * This probably only works for water!!! PvM
     */
    if ((xcm[XX] - TOLERANCE) > 0.0) {
      make_rot_mat(ZZ,-0.5*M_PI,r_mat);
      rotate_mol(k0,k1,index,x,r_mat);
    }
  }
}

/* Does the real work */
void do_multipoles(char *trjfn,char *topfn,char *molndxfn,gmx_bool bFull)
{
  int        i;
  int        gnx;
  atom_id    *grpindex;
  char       *grpname;
  t_topology *top;
  t_atom     *atom;
  t_block    *mols;
  int        natoms,status;
  real       t;
  matrix     box;
  real       t0,t1,tq;
  int        teller;
  gmx_bool       bCont;

  rvec       *x,*m1;
  tensor     *m2;
  tensor3    *m3;
  tensor4    *m4;
  matrix     trans;
  gmx_rmpbc_t  gpbc=NULL;

  top  = read_top(topfn);
  rd_index(molndxfn,1,&gnx,&grpindex,&grpname);
  atom = top->atoms.atom;
  mols = &(top->blocks[ebMOLS]);

  natoms  = read_first_x(&status,trjfn,&t,&x,box);
  snew(m1,gnx);
  snew(m2,gnx);
  snew(m3,gnx);
  snew(m4,gnx);

  gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);

  /* Start while loop over frames */
  do {
    /* PvM, bug in rm_pbc??? Does not work for virtual sites ...
    gmx_rmpbc(gpbc,box,x,x_s); */

    /* Begin loop of all molecules in index file */
    for(i=0; (i<gnx); i++) {
      int gi = grpindex[i];

      rot_mol_to_std_orient(mols->index[gi],mols->index[gi+1],mols->a,atom,x,
			    trans);

      /* Rotate the molecule along the principal moments axes */
      /* rotate_mol(mols->index[gi],mols->index[gi+1],mols->a,x,trans); */

      /* Compute the multipole moments */
      mol_M1(mols->index[gi],mols->index[gi+1],mols->a,x,atom,m1[i]);
      mol_M2(mols->index[gi],mols->index[gi+1],mols->a,x,atom,m2[i]);
      mol_M3(mols->index[gi],mols->index[gi+1],mols->a,x,atom,m3[i]);
      mol_M4(mols->index[gi],mols->index[gi+1],mols->a,x,atom,m4[i]);

      /* Spit it out */
      pr_M1(stdout,"",i,m1[i],t);
      pr_M2(stdout,"",m2[i],bFull);
      pr_M3(stdout,"",m3[i],bFull);
      pr_M4(stdout,"",m4[i],bFull);


    } /* End loop of all molecules in index file */
    
    bCont = read_next_x(status,&t,natoms,x,box);
  } while(bCont);
  gmx_rmpbc_done(gpbc);

  
  
}

int gmx_multipoles(int argc,char *argv[])
{
  const char *desc[] = {
    "g_multipoles computes the electric multipole moments of",
    "molecules selected by a molecular index file.",
    "The center of mass of the molecule is used as the origin"
  };

  static gmx_bool bFull = FALSE;
  static int  ntb=0;
  t_pargs pa[] = {
    { "-boxtype",FALSE,etINT,&ntb, "HIDDENbox type 0=rectangular; 1=truncated octahedron (only rectangular boxes are fully implemented)"},
    { "-full",   FALSE, etBOOL, &bFull, 
      "Print all compononents of all multipoles instead of just the interesting ones" }
  };
  int          gnx;
  atom_id      *index;
  char         *grpname;
  t_filenm fnm[] = {
    { efTPX, NULL, NULL,      ffREAD },
    { efTRX, "-f", NULL,      ffREAD },
    { efNDX, NULL, NULL,      ffREAD }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  int         i,j,k;
  int         natoms;
  int         step;
  real        t,lambda;

  t_tpxheader tpx;
  t_topology  top;
  rvec        *x,*xnew;
  matrix      box;
  t_atom      *atom;
  rvec        dipole,dipole2;
  real        mtot,alfa,beta,gamma;
  rvec        CoM,*xCoM,angle;
  real        *xSqrCoM;

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  do_multipoles(ftp2fn(efTRX,NFILE,fnm),ftp2fn(efTPX,NFILE,fnm),
		ftp2fn(efNDX,NFILE,fnm),bFull);

  return 0;
}
