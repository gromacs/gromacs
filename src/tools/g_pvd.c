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
 * Grunge ROck MAChoS
 */
#include <ctype.h>
#include "sysstuff.h"
#include "statusio.h"
#include "vec.h"
#include "smalloc.h"
#include "physics.h"
#include "xvgr.h"
#include "futil.h"
#include "maths.h"
#include "rdgroup.h"
#include "statutil.h"
#include "copyrite.h"

#define MVBUF  1000000

void plot_vol(FILE *status, char *volout)
{
  float        volbuf[MVBUF];
  FILE         *out;
  t_statheader sh;
  int          ct,idum,i,j,nframes=0;
  real         rdum,t,vol,voltot=0;
  double       vtmp,vav2,vav,vrms;
  matrix       box;

  out=xvgropen(volout,"Volume of Computational Box",
	       "Time (ps)","Volume (nm\\S3\\N)");
  j=0;
  while (!eof(status)) {
    rd_header(status,&sh);
    rd_hstatus(status,&sh,&idum,&t,&rdum,NULL,
	       sh.box_size ? box : NULL,NULL,NULL,&idum,NULL,NULL,NULL,
	       &idum,NULL,NULL);
    
    ct=check_times(t);
    
    if (ct > 0)
      break;
    else if (ct == 0) {
      if (sh.box_size) {
        vol=det(box);
        voltot+=vol;
        volbuf[nframes]=vol;
        nframes++;
        fprintf(out,"%12.5e  %12.5e\n",t,vol);
        if ((j++ % 10) == 0)
	  fprintf(stderr,"\rFrame: %5d",j-1);
      }
    }
  }
  fprintf(stderr,"\n");
 
  vav=voltot/nframes;
  vav2=0;
  for(i=0; (i<nframes); i++) {
    vtmp=(volbuf[i]-vav);
    vav2+=vtmp*vtmp;
  }
  vrms=sqrt(vav2/nframes);
  fprintf(stderr,"Average Volume: %10g RMS: %10g\n",vav,vrms);
  fprintf(out,"@    subtitle \"Average: %g, RMS: %g\"\n",vav,vrms);
  fclose(out);
  
}

void plot_box(FILE *status, char *volout)
{
  static char *nms[] = { "X", "Y", "Z" };
  FILE         *out;
  t_statheader sh;
  int          j,idum,ct;
  real         t,rdum;
  matrix       box;

  out=xvgropen(volout,"Box Length",
	       "Time (ps)","Length (nm)");
  xvgr_legend(out,asize(nms),nms);
  j=0;
  while (!eof(status)) {
    rd_header(status,&sh);
    rd_hstatus(status,&sh,&idum,&t,&rdum,NULL,
	       sh.box_size ? box : NULL,NULL,NULL,&idum,NULL,NULL,NULL,
	       &idum,NULL,NULL);
    
    ct=check_times(t);
    
    if (ct > 0)
      break;
    else if (ct == 0) {
      if (sh.box_size) {
        fprintf(out,"%12.5e  %12.5e  %12.5e  %12.5e\n",
		t,box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
        if ((j++ % 10) == 0)
	  fprintf(stderr,"\rFrame: %5d",j-1);
      }
    }
  }
  fprintf(stderr,"\n");
 
  fclose(out);
}

void plot_pres(FILE *status, char *presout)
{
  static char *ptitle[] = { "P-x", "P-y", "P-z", "P-total" };
  FILE         *out;
  t_statheader sh;
  int          idum,j;
  real         rdum,t;
  tensor       pres;
  real         tp;
  int          np,ct;

  out=xvgropen(presout,"Gromacs Pressures","Time (ps)","P (bar)");
  xvgr_legend(out,asize(ptitle),ptitle);
  j=np=0;
  while (!eof(status)) {
    rd_header(status,&sh);
    rd_hstatus(status,&sh,&idum,&t,&rdum,NULL,
	       NULL,NULL,sh.pres_size ? pres : NULL,&idum,NULL,NULL,NULL,
	       &idum,NULL,NULL);
    ct=check_times(t);
    if (ct > 0)
      break;
    else if ((ct == 0) && (sh.pres_size > 0)) {
      if ((j++ % 10) == 0)
	fprintf(stderr,"\rFrame: %5d",j);
      tp=trace(pres);
      fprintf(out,"%10g  %10g  %10g  %10g  %10g\n",
	      t,
	      PRESFAC*pres[XX][XX],
	      PRESFAC*pres[YY][YY],
	      PRESFAC*pres[ZZ][ZZ],
	      PRESFAC*tp/3.0);
    }
  }
  fclose(out);
}

void plot_prestens(FILE *status, char *presout, char *Qname)
{

  FILE         *out,*out2;
  t_statheader sh;
  int          idum,i,j,d,k,fr;
  real         rdum,t,t_0,var_P;
  tensor       pres,pres_a;

  if ((out=fopen(presout,"w"))==NULL) {
    perror(presout);
    exit(1);
  }
  out2=fopen(Qname,"w");

  rd_header(status,&sh);
  rd_hstatus(status,&sh,&idum,&t,&rdum,NULL,
             NULL,NULL,sh.pres_size ? pres : NULL,&idum,NULL,NULL,NULL,
             &idum,NULL,NULL);
  t_0=t;
  rewind(status);

  fr=0;
  while (!eof(status)) {
    rd_header(status,&sh);
    rd_hstatus(status,&sh,&idum,&t,&rdum,NULL,
	       NULL,NULL,sh.pres_size ? pres : NULL,&idum,NULL,NULL,NULL,
	       &idum,NULL,NULL);
    fr++;
    if (sh.pres_size) {
      fprintf(stderr,"\rReading Frame: %5d",fr);
      fflush(stderr);
      fprintf(out,"%12.5e  ",t);
      for (d=0;(d<DIM); d++){
        for (k=0;(k<DIM); k++){
          fprintf(out,"%12.5e  ", pres[d][k]);
        }
      }
      fprintf(out,"\n");


/*      fprintf(out2,"%12.5e  %12.5e  %12.5e  %12.5e  \n",t,pres[XX][YY]/pres[YY][XX],pres[XX][ZZ]/pres[ZZ][XX],pres[YY][ZZ]/pres[ZZ][YY]); */

      
      /* calculate the anti-symmetric pressure tensor elements */
      var_P=0;
      for (i=0; (i<DIM); i++){
        for (j=0; (j<DIM); j++){
          pres_a[i][j]= pres[i][j]-pres[j][i];
          var_P += pres_a[i][j]*pres[i][j];
        }
      }
      /* normalisation */
      var_P /= DIM*(DIM-1);
      fprintf(out2,"%12.5e  %12.5e\n",t,sqrt(var_P));
    }
  }
  /* xvgr stuff */
  fprintf(out2,"@ yaxis label \"Root mean square P\\SA\\N elements (bar)\"\n");
  fprintf(out2,"@ xaxis label \"Time (ps)\"\n");
  fprintf(out2,"@ title \"%s\"\n",status_title(status));


  fprintf(out,"time  Pxx  Pxy  Pxz  Pyx  Pyy  Pyz  Pzx  Pzy  Pzz"); 
  fprintf(out,"TitleText: Gromacs: %s\n",status_title(status));
  fprintf(out,"XUnitText: Time (ps)\n");
  fprintf(out,"YUnitText: Pressure (bar)\n");
  fclose(out);
  fclose(out2);
  fprintf(stderr,"\n");
}

void plot_density(FILE *status,char *densout,int axis,
		  int ngrp,
		  int nx[],atom_id *index[],char *gname[])
{
#define NPX 51
  typedef int ivec[DIM];
  static char *szAxis[DIM] = { "X (nm)", "Y (nm)", "Z (nm)" };
  FILE   *out;
  real   t;
  int    i,ii,ct,ind,j,natoms;
  rvec   *x;
  matrix box;
  real   nframes,slab;
  real   seg,invseg;
  int    iseg;
  real   **px;

  out=xvgropen(densout,"Gromacs Density Profiles",szAxis[axis],
	       "Density (atoms/nm\\S3\\N)");
  xvgr_legend(out,ngrp,gname);
  
  if ((natoms=read_first_x(status,&t,&x,box))==0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }
  /* Rectangular box only... */
  snew(px,ngrp);
  for(ii=0; (ii<ngrp); ii++)
    snew(px[ii],NPX);
  seg=box[axis][axis]/(NPX-1.0);
  invseg=1.0/seg;
  j=0;
  do {
    ct=check_times(t);
    if (ct > 0)
      break;
    else if (ct == 0) {
      if ((j++ % 10) == 0)
	fprintf(stderr,"\rframe: %5d",j-1);
      
      /* Calculate box-segment */
      for(ii=0; (ii<ngrp); ii++) {
	for(i=0; (i<nx[ii]); i++) {
	  ind=index[ii][i];
	  iseg=gmx_nint(x[ind][axis]*invseg);
	  px[ii][iseg]+=1;
	}
      }
    }
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  /* Calculate the volume per box-slab */
  slab=det(box);
  slab/=(NPX-1.0);

  nframes=1.0/j;
  for(ii=0; (ii<ngrp); ii++) {
    px[ii][0]+=px[ii][NPX-1];
  }
  for(i=0; (i<NPX-1); i++) {
    fprintf(out,"%10g  ",i*seg);
    for(ii=0; (ii<ngrp); ii++) {
      px[ii][i]*=nframes;
      fprintf(out,"%10g  ",px[ii][i]);
    }
    fprintf(out,"\n");
  }

  fclose(out);
  
}

enum { eBox, eVol, ePres, eDens ,ePresTens };

static int get_funct(int argc,char *argv[],int *axis)
{
  int ret=0;

  if (argc < 2)
    usage(argv[0],argv[0]);
  switch (argv[1][1]) {
  case 'x':
    ret=eBox;
    break;
  case 'v':
    ret=eVol;
    break;
  case 'p':
    ret=ePres;
    break;
  case 'd':
    ret=eDens;
    if (argc < 3)
      usage(argv[0],argv[1]);
    *axis=toupper(argv[2][0])-'X';
    break;
  default:
    usage(argv[0],argv[0]);
  }
  return ret;
}

void main(int argc,char *argv[])
{
  static char *desc[] = {
    "Read a trajectory and plot either volume, pressure, density",
    "or box edges as a function of time."
  };
  static char *opts[] = {
    " -f",
    " -x",
    " -v",
    " -p",
    " -P",
    " -d X | Y | Z"
  };
  static char *odesc[] = {
    "trajectory",
    "plot box lengths",
    "plot volume",
    "plot pressure (px,py,pz,p)",
    "plot complete pressure tensor P ; write root mean square of off-diagonalelements",
    "plot density in X, Y or Z direction of various components of the system using an indexfile"
  };
  t_manual man = {asize(desc),desc,asize(opts),opts,odesc,0,NULL};
  FILE      *status;
  int       funct;
  char      **grpname;
  int       ngrp,*gnx,axis;
  atom_id   **index;
  t_filenm fnm[] = {
    { efTRJ, "-f", NULL, ffREAD },
    { efNDX, NULL, NULL, ffREAD },
    { efXVG, NULL, NULL, ffWRITE },
    { efXVG, "-P", "Q",  ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,
		    NFILE,fnm,TRUE,&man);
  funct=get_funct(argc,argv,&axis);
  
  status=ftp2FILE(efTRJ,NFILE,fnm,"r");
  switch (funct) {
  case eBox:
    plot_box(status,ftp2fn(efXVG,NFILE,fnm));
    break;
  case eVol:
    plot_vol(status,ftp2fn(efXVG,NFILE,fnm));
    break;
  case ePres:
    plot_pres(status,ftp2fn(efXVG,NFILE,fnm));
    break;
  case ePresTens:
    plot_prestens(status,ftp2fn(efXVG,NFILE,fnm),opt2fn("-P",NFILE,fnm));
    break;
  case eDens:
    printf("How many groups ? ");
    do { scanf("%d",&ngrp); } while (ngrp <= 0);
    snew(grpname,ngrp);
    snew(index,ngrp);
    snew(gnx,ngrp);
    rd_index(ftp2fn(efNDX,NFILE,fnm),ngrp,gnx,index,grpname);
    plot_density(status,ftp2fn(efXVG,NFILE,fnm),axis,
		 ngrp,gnx,index,grpname);
    break;
  }
  fclose(status);
  
  if (funct != ePresTens)
    xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
  thanx(stdout);
}
