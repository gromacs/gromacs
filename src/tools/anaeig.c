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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#include <math.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "confio.h"
#include "physics.h"
#include "readev.h"

typedef struct {
  real x,y,z;
  int  i;
} t_point;

int ptcomp(const void *a,const void *b)
{
  t_point *av,*bv;
  
  av=(t_point *)a;
  bv=(t_point *)b;
  
  if (av->z < bv->z)
    return -1;
  else if (av->z > bv->z)
    return 1;
  return 0;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "anaeig analyzes eigenvector from essential dynamics."
  };
  static char *opts[] = {
    "-ev"
  };
  static char *odesc[] = {
    "Eigenvector index (default 1)"
  };
  t_manual man = { asize(desc),desc,asize(opts),opts,odesc,0,NULL};
  t_filenm  fnm[] = {
    { efGRO, "-c", NULL,      FALSE },
    { efDAT, "-d","eigenvec", FALSE },
    { efXVG, "-o","decomp",   FALSE },
    { efXVG, "-p","plane",    FALSE }
  };
#define NFILE asize(fnm)
  static char *setname[] = {
    "X", "Y", "Z (Axial)", "Radial", "Tangential"
  };
  
  FILE    *out;
  char    buf[256];
  rvec    *x,*v;
  t_point *pt;
  rvec    **EV;
  real    aa,alpha,ca,sa,*evx,*evy;
  matrix  box;
  int     natoms,ev=1;
  int     i,j,ei,ai;
  double  xx,yy,zz,rad,phi0,phi1,dphi,r0,r1,r2,dr,x0,y0,x1,y1;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,
		    NFILE,fnm,TRUE,&man);
  
  for(i=1; (i<argc); i++) {
    if (i==argc-1)
      usage(argv[0],argv[i]);
    else if (strcmp(argv[i],"-ev") == 0) {
      ev=atoi(argv[++i]);
    }
  }
  fprintf(stderr,"Will decompose eigenvector %d\n",ev);
  
  get_coordnum(ftp2fn(efGRO,NFILE,fnm),&natoms);
  snew(x,natoms);
  snew(v,natoms);
  read_conf(ftp2fn(efGRO,NFILE,fnm),buf,&natoms,x,v,box);

  EV=read_ev(ftp2fn(efDAT,NFILE,fnm),natoms);
  
  sprintf(buf,"Decomposition of EV %d",ev);
  out=xvgropen(opt2fn("-o",NFILE,fnm),buf,"Atom","nm");
  xvgr_legend(out,asize(setname),setname);  

  for(j=0; (j<natoms); j++) {
    xx=EV[ev-1][j][XX];
    yy=EV[ev-1][j][YY];
    zz=EV[ev-1][j][ZZ];
    
    x0=x[j][XX];
    y0=x[j][YY];
    x1=x0+xx;
    y1=y0+yy;
    r0=sqrt(sqr(x0)+sqr(y0));
    r1=sqrt(sqr(x1)+sqr(y1));
    rad=r0-r1;
    
    phi0=atan2(y0,x0);
    phi1=atan2(y1,x1);
    dphi=RAD2DEG*(phi1-phi0)+360;
    while (dphi > 180)
      dphi-=360;
    while (dphi < -180)
      dphi+=360;
    
    r2=sqr(xx)+sqr(yy);
    dr=sqrt(r2-sqr(rad));
    if (dphi < 0)
      dr*=-1;
    
    fprintf(out,"%5d  %10g  %10g  %10g  %10g  %10g\n",
	    j,xx,yy,zz,rad,dr);
  }
  fclose(out);

  /* Now rotate eigenvector for bending test */
  snew(pt,natoms);
  snew(evx,natoms);
  snew(evy,natoms);
  for(i=0; (i<natoms); i++) {
    evx[i]=EV[ev-1][i][XX];
    evy[i]=EV[ev-1][i][YY];
  }
  lsq_y_ax(natoms,evx,evy,&aa);
  alpha=atan(aa);
  ca=cos(alpha);
  sa=sin(alpha);
  for(i=0; (i<natoms); i++) {
    pt[i].x= ca*evx[i]+sa*evy[i];
    pt[i].y=-sa*evx[i]+ca*evy[i];
    pt[i].z=EV[ev-1][i][ZZ];         /*for sorting: x[i][ZZ]*/
    pt[i].i=i;
  }
  /*qsort((char *)pt,natoms,sizeof(pt[0]),ptcomp);*/
  sprintf(buf,"EV %d rotated %.0f deg around Z-axis",ev,RAD2DEG*alpha);
  out=xvgropen(opt2fn("-p",NFILE,fnm),buf,"Atom","nm");
  for(i=0; (i<natoms); i++) {
    fprintf(out,"%10d  %10g  %10g  %10g\n",pt[i].i,pt[i].x,pt[i].y,pt[i].z);
  }
  fclose(out);

  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
  		      
  thanx(stdout);
  
  return 0;
}

