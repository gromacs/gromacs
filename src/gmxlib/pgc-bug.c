/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GROup of MAchos and Cynical Suckers
 */
#include <stdio.h>

#define XX  0
#define YY  1
#define ZZ  2
#define DIM 3
typedef float    real;
typedef real     rvec [DIM];
typedef rvec     tensor[DIM];
typedef unsigned int  t_ishift;
	
static void upd_vir(rvec vir,real dvx,real dvy,real dvz)
{
  vir[XX]-=0.5*dvx;
  vir[YY]-=0.5*dvy;
  vir[ZZ]-=0.5*dvz;
}

static void lo_fcv(int i0,int i1,int g0,
		   real x[],real f[],tensor vir,
		   t_ishift is[],real shift_vec[])
{
  int      i,i3,gg,t,t3;
  real     xx,yy,zz;
  real     dvxx=0,dvxy=0,dvxz=0,dvyx=0,dvyy=0,dvyz=0,dvzx=0,dvzy=0,dvzz=0;

  for(i=i0,gg=g0; (i<i1); i++,gg++) {
    i3=DIM*i;
    t=is[gg];
    t3=DIM*t;
    
    xx=x[i3+XX]-shift_vec[t3+XX];
    dvxx+=xx*f[i3+XX];
    dvxy+=xx*f[i3+YY];
    dvxz+=xx*f[i3+ZZ];
    
    yy=x[i3+YY]-shift_vec[t3+YY];
    dvyx+=yy*f[i3+XX];
    dvyy+=yy*f[i3+YY];
    dvyz+=yy*f[i3+ZZ];
    
    zz=x[i3+ZZ]-shift_vec[t3+ZZ];
    dvzx+=zz*f[i3+XX];
    dvzy+=zz*f[i3+YY];
    dvzz+=zz*f[i3+ZZ];
  }
  
  upd_vir(vir[XX],dvxx,dvxy,dvxz);
  upd_vir(vir[YY],dvyx,dvyy,dvyz);
  upd_vir(vir[ZZ],dvzx,dvzy,dvzz);
}

int main(int argc,char *argv[])
{
#define MAX 100
  int i0=10,i1=50;
  int g0=10;
  rvec x[MAX];
  rvec f[MAX];
  tensor vir;
  t_ishift is[MAX];
   rvec shift_vec[MAX];
  int i;

  for(i=0; (i<MAX); i++) {
    x[i][XX]=x[i][YY]=x[i][ZZ]=0;
    f[i][XX]=f[i][YY]=f[i][ZZ]=0;
    shift_vec[i][XX]=shift_vec[i][YY]=shift_vec[i][ZZ]=0;
    is[i]=0;
  }
  for(i=0; (i<DIM); i++)
    vir[i][XX]=vir[i][YY]=vir[i][ZZ]=0;
    
  lo_fcv(i0,i1,g0,x[0],f[0],vir,is,shift_vec[0]);
}
