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
static char *SRCID_g_mindist_c = "$Id$";

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
#include "rdgroup.h"

static void calc_mindist(real MinDist,
			 matrix box, rvec x[], 
			 int nx1,int nx2,
			 atom_id index1[], atom_id index2[],
			 real *md, int *nd,
			 int *ixmin, int *jxmin)
{
  int     i,j,j0=0,j1;
  int     ix,jx;
  atom_id *index3;
  rvec    dx;
  int     numd=0;
  real    r2,rmin2,r_2;

  *ixmin = -1;
  *jxmin = -1;
  
/*   real    hb2; */
/*   hb2=(sqr(box[XX][XX])+sqr(box[YY][YY])+sqr(box[ZZ][ZZ]))*0.5; */
  
  r_2=sqr(MinDist);

  /* Must init pbc every step because of pressure coupling */
  init_pbc(box,FALSE);
  if (index2) {
    j0=0;
    j1=nx2;
    index3=index2;
  }
  else {
    j1=nx1;
    index3=index1;
  }

  rmin2=1000.0;

  for(i=0; (i < nx1); i++) {
    ix=index1[i];
    if (!index2)
      j0=i+1;
    for(j=j0; (j < j1); j++) {
      jx=index3[j];
      if (ix != jx) {
	pbc_dx(x[ix],x[jx],dx);
	r2=iprod(dx,dx);
	if (r2 < rmin2)
	  rmin2=r2;
	if (r2 < r_2) {
	  numd++;
	  *ixmin=ix;
	  *jxmin=jx;
	}
      }
    }
  }
  *md=sqrt(rmin2);
  *nd=numd;
}

void mindist_plot(char *fn,FILE *atm,real mind,
		  char *dfile,char *nfile,bool bMat,
		  int ng,atom_id *index[], int gnx[], char *grpn[])
{
  FILE         *dist,*num;
  char         buf[256];
  char         **leg;
  real         t,md;
  int          nd,status;
  int          i,j,k,natoms;
  int	       min1,min2;
  rvec         *x0;
  matrix       box;
  
  if ((natoms=read_first_x(&status,fn,&t,&x0,box))==0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }

  sprintf(buf,"Number of Contacts < %g nm",mind);
  dist=xvgropen(dfile,"Minimum Distance","Time (ps)","Distance (nm)");
  num=xvgropen(nfile,buf,"Time (ps)","Number");

  if (bMat) {
    snew(leg,(ng*(ng-1))/2);
    for(i=j=0; (i<ng-1); i++) {
      for(k=i+1; (k<ng); k++,j++) {
	sprintf(buf,"%s-%s",grpn[i],grpn[k]);
	leg[j]=strdup(buf);
      }
    }
    xvgr_legend(dist,j,leg);
    xvgr_legend(num,j,leg);
  }
  else {  
    snew(leg,ng-1);
    for(i=0; (i<ng-1); i++) {
      sprintf(buf,"%s-%s",grpn[0],grpn[i+1]);
      leg[i]=strdup(buf);
    }
    xvgr_legend(dist,ng-1,leg);
    xvgr_legend(num,ng-1,leg);
  }    
  j=0;
  do {
    if ((j++ % 10) == 0)
      fprintf(stderr,"\rframe: %5d",j-1);
      
    fprintf(dist,"%12g",t);
    fprintf(num,"%12g",t);

    if (bMat) {
      for(i=0; (i<ng-1); i++) {
	for(k=i+1; (k<ng); k++) {
	  calc_mindist(mind,box,x0,gnx[i],gnx[k],
		       index[i],index[k],&md,&nd,
		       &min1,&min2);
	  fprintf(dist,"  %12g",md);
	  fprintf(num,"  %8d",nd);
	}
      }
    }
    else {    
      for(i=1; (i<ng); i++) {
	calc_mindist(mind,box,x0,gnx[0],gnx[i],
		     index[0],index[i],&md,&nd,
		     &min1,&min2);
	fprintf(dist,"  %12g",md);
	fprintf(num,"  %8d",nd);
      }    
    }
    fprintf(dist,"\n");
    fprintf(num,"\n");
    if (min1 != -1)
      fprintf(atm,"%12g  %12d  %12d\n",t,min1,min2);
  } while (read_next_x(status,&t,natoms,x0,box));
  
  fprintf(stderr,"\n");
  close_trj(status);
  fclose(dist);
  fclose(num);

  sfree(x0);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_mindist computes the distance between one group and a number of",
    "other groups.",
    "Both the smallest distance and the number of contacts within a given",
    "distance are plotted to two separate output files"
  };
  
  static bool bMat=FALSE;
  static real mindist=0.6;
  t_pargs pa[] = {
    { "-matrix", FALSE, etBOOL, &bMat,
      "Calculate half a matrix of group-group distances" },
    { "-d",      FALSE, etREAL, &mindist,
      "Distance for contacts" }
  };
  FILE      *status, *atm;
  int       i,ng;
  char      **grpname;
  int       *gnx;
  atom_id   **index;
  t_filenm  fnm[] = {
    { efTRX, "-f", NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffREAD },
    { efXVG, "-od","dist",ffWRITE },
    { efXVG, "-on","num", ffWRITE },
    { efOUT, "-o","atm-pair", ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
		    
  if (bMat)
    fprintf(stderr,"You can compute all distances between a number of groups\n"
	    "How many groups do you want (>= 2) ?\n");
  else
    fprintf(stderr,"You can compute the distances between a first group\n"
	    "and a number of other groups.\n"
	    "How many other groups do you want (>= 1) ?\n");
  ng=0;
  do {
    scanf("%d",&ng);
    if (!bMat)
      ng++;
  } while (ng < 2);
  snew(gnx,ng);
  snew(index,ng);
  snew(grpname,ng);
  
  rd_index(ftp2fn(efNDX,NFILE,fnm),ng,gnx,index,grpname);

  atm=ftp2FILE(efOUT,NFILE,fnm,"w");
  mindist_plot(ftp2fn(efTRX,NFILE,fnm),atm,mindist,
	       opt2fn("-od",NFILE,fnm),opt2fn("-on",NFILE,fnm),
	       bMat,ng,index,gnx,grpname);
  fclose(atm);

  xvgr_file(opt2fn("-od",NFILE,fnm),"-nxy");
  xvgr_file(opt2fn("-on",NFILE,fnm),"-nxy");
  
  thanx(stdout);
  
  return 0;
}

