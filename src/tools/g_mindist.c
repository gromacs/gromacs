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
 * GRowing Old MAkes el Chrono Sweat
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
#include "tpxio.h"
#include "rmpbc.h"

static void periodic_dist(matrix box,rvec x[],int n,atom_id index[],
			  real *rmin,real *rmax)
{
#define NSHIFT 26
  int  sx,sy,sz,i,j,s;
  real sqr_box,r2min,r2max,r2;
  rvec shift[NSHIFT],d0,d;

  sqr_box = sqr(min(box[XX][XX],min(box[YY][YY],box[ZZ][ZZ])));

  s = 0;
  for(sz=-1; sz<=1; sz++)
    for(sy=-1; sy<=1; sy++)
      for(sx=-1; sx<=1; sx++)
	if (sx!=0 || sy!=0 || sz!=0) {
	  for(i=0; i<DIM; i++)
	    shift[s][i] = sx*box[XX][i]+sy*box[YY][i]+sz*box[ZZ][i];
	  s++;
	}
  
  r2min = sqr_box;
  r2max = 0;

  for(i=0; i<n; i++)
    for(j=i+1; j<n; j++) {
      rvec_sub(x[index[i]],x[index[j]],d0);
      r2 = norm2(d0);
      if (r2 > r2max)
	r2max = r2;
      for(s=0; s<NSHIFT; s++) {
	rvec_add(d0,shift[s],d);
	r2 = norm2(d);
	if (r2 < r2min)
	  r2min = r2;
      }
    }

  *rmin = sqrt(r2min);
  *rmax = sqrt(r2max);
}

static void periodic_mindist_plot(char *trxfn,char *outfn,
				  t_topology *top,int n,atom_id index[])
{
  FILE   *out;
  char   *leg[5] = { "min per.","max int.","box1","box2","box3" };
  int    status;
  real   t;
  rvec   *x;
  matrix box;
  int    natoms;
  real   r,rmin,rmax,rmint,tmint;

  natoms=read_first_x(&status,trxfn,&t,&x,box);
  
  check_index(NULL,n,index,NULL,natoms);

  out = xvgropen(outfn,"Minimum distance to periodic image",
		 "Time (ps)","Distance (nm)");
  fprintf(out,"@ subtitle \"and maximum internal distance\"\n");
  xvgr_legend(out,5,leg);

  rmint = box[XX][XX];
  tmint = 0;

  do {
    rm_pbc(&(top->idef),natoms,box,x,x);
    periodic_dist(box,x,n,index,&rmin,&rmax);
    if (rmin < rmint) {
      rmint = rmin;
      tmint = t;
    }
    fprintf(out,"\t%g\t%6.3f %6.3f %6.3f %6.3f %6.3f\n",t,rmin,rmax,
	    norm(box[0]),norm(box[1]),norm(box[2]));
  } while(read_next_x(status,&t,natoms,x,box));

  fclose(out);

  fprintf(stdout,
	  "\nThe shortest periodic distance is %g (nm) at time %g (ps)\n",
	  rmint,tmint);
}

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
  
  if ((natoms=read_first_x(&status,fn,&t,&x0,box))==0)
    fatal_error(0,"Could not read coordinates from statusfile\n");

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
    "Both the minimum distance and the number of contacts within a given",
    "distance are written to two separate output files.[PAR]",
    "With option [TT]-pi[tt] the minimum distance of a group to its",
    "periodic image is plotted. This is useful for checking if a protein",
    "has seen its periodic image during a simulation. Only one shift in",
    "each direction is considered, giving a total of 26 shifts.",
    "It also plots the maximum distance within the group and the lengths",
    "of the three box vectors.",
    "This option is very slow."
  };
  
  static bool bMat=FALSE,bPer=FALSE;
  static real mindist=0.6;
  t_pargs pa[] = {
    { "-matrix", FALSE, etBOOL, {&bMat},
      "Calculate half a matrix of group-group distances" },
    { "-d",      FALSE, etREAL, {&mindist},
      "Distance for contacts" },
    { "-pi",      FALSE, etBOOL, {&bPer},
      "Calculate minimum distance with periodic images" }
  };
  t_topology top;
  char       title[256];
  real       t;
  rvec       *x;
  matrix     box;
  
  FILE      *atm;
  int       ng;
  char      *tps,*ndx,**grpname;
  int       *gnx;
  atom_id   **index;
  t_filenm  fnm[] = {
    { efTRX, "-f", NULL,  ffREAD },
    { efTPS, NULL, NULL,  ffOPTRD },
    { efNDX, NULL, NULL,  ffOPTRD },
    { efXVG, "-od","mindist",ffWRITE },
    { efXVG, "-on","numcont", ffWRITE },
    { efOUT, "-o","atm-pair", ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  tps = ftp2fn_null(efTPS,NFILE,fnm);
  ndx = ftp2fn_null(efNDX,NFILE,fnm);
  
  if (bPer) {
    ng = 1;
    fprintf(stderr,"Choose a group for distance calculation\n");
  } else {
    if (bMat)
      fprintf(stderr,"You can compute all distances between a number of groups\n"
	      "How many groups do you want (>= 2) ?\n");
    else
      fprintf(stderr,"You can compute the distances between a first group\n"
	      "and a number of other groups.\n"
	      "How many other groups do you want (>= 1) ?\n");
    ng = 0;
    do {
      scanf("%d",&ng);
      if (!bMat)
	ng++;
    } while (ng < 2);
  }
  snew(gnx,ng);
  snew(index,ng);
  snew(grpname,ng);

  if (tps || ndx==NULL)
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&x,NULL,box,FALSE);
  
  get_index(&top.atoms,ndx,ng,gnx,index,grpname);
  
  if (bPer) {
    periodic_mindist_plot(ftp2fn(efTRX,NFILE,fnm),opt2fn("-od",NFILE,fnm),
			  &top,gnx[0],index[0]);
  } else {
    atm=ftp2FILE(efOUT,NFILE,fnm,"w");
    mindist_plot(ftp2fn(efTRX,NFILE,fnm),atm,mindist,
		 opt2fn("-od",NFILE,fnm),opt2fn("-on",NFILE,fnm),
		 bMat,ng,index,gnx,grpname);
    fclose(atm);
  }

  xvgr_file(opt2fn("-od",NFILE,fnm),"-nxy");
  if (!bPer)
    xvgr_file(opt2fn("-on",NFILE,fnm),"-nxy");
  
  thanx(stderr);
  
  return 0;
}

