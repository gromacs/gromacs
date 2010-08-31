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
#include <string.h>

#include "macros.h"
#include "vec.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "filenm.h"
#include "statutil.h"
#include "copyrite.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "pbc.h"
#include "xvgr.h"
#include "gmx_ana.h"


typedef struct {
  char *label;
  int  cg;
  real q;
} t_charge;

static t_charge *mk_charge(t_atoms *atoms,t_block *cgs,int *nncg)
{
  t_charge *cg=NULL;
  char     buf[32];
  int      i,j,ncg,resnr,anr;
  real     qq;

  /* Find the charged groups */  
  ncg=0;
  for(i=0; (i<cgs->nr); i++) {
    qq=0.0;
    for(j=cgs->index[i]; (j<cgs->index[i+1]); j++) {
      qq+=atoms->atom[j].q;
    }
    if (fabs(qq) > 1.0e-5) {
      srenew(cg,ncg+1);
      cg[ncg].q=qq;
      cg[ncg].cg=i;
      anr=cgs->index[i];
      resnr=atoms->atom[anr].resind;
      sprintf(buf,"%s%d-%d",
	      *(atoms->resinfo[resnr].name),
	      atoms->resinfo[resnr].nr,
	      anr+1);
      cg[ncg].label=strdup(buf);
      ncg++;
    }
  }
  *nncg=ncg;
  
  for(i=0; (i<ncg); i++) {
    printf("CG: %10s Q: %6g  Atoms:",
	   cg[i].label,cg[i].q);
    for(j=cgs->index[cg[i].cg]; (j<cgs->index[cg[i].cg+1]); j++)
      printf(" %4u",j);
    printf("\n");
  }
  
  return cg;
}

static real calc_dist(t_pbc *pbc,rvec x[],t_block *cgs,int icg,int jcg)
{
  int  i,j;
  rvec dx;
  real d2,mindist2=1000;

  for(i=cgs->index[icg]; (i<cgs->index[icg+1]); i++) {
    for(j=cgs->index[jcg]; (j<cgs->index[jcg+1]); j++) {
      pbc_dx(pbc,x[i],x[j],dx);
      d2 = norm2(dx);
      if (d2 < mindist2)
	mindist2 = d2;
    }
  }
  return sqrt(mindist2);
}

int gmx_saltbr(int argc,char *argv[])
{
  const char *desc[] = {
    "g_saltbr plots the distance between all combination of charged groups",
    "as a function of time. The groups are combined in different ways.",
    "A minimum distance can be given, (eg. the cut-off), then groups",
    "that are never closer than that distance will not be plotted.[BR]",
    "Output will be in a number of fixed filenames, min-min.xvg, plus-min.xvg",
    "and plus-plus.xvg, or files for every individual ion-pair if the [TT]-sep[tt]",
    "option is selected. In this case files are named as [TT]sb-ResnameResnr-Atomnr[tt].",
    "There may be many such files."
  };
  static gmx_bool bSep=FALSE;
  static real truncate=1000.0;
  t_pargs pa[] = {
    { "-t",   FALSE, etREAL, {&truncate},
      "trunc distance" },
    { "-sep", FALSE, etBOOL, {&bSep},
      "Use separate files for each interaction (may be MANY)" }
  };
  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL, ffREAD },
    { efTPX, NULL,  NULL, ffREAD },
  };
#define NFILE asize(fnm)

  FILE       *out[3],*fp;
  static const char *title[3] = {
    "Distance between positively charged groups",
    "Distance between negatively charged groups",
    "Distance between oppositely charged groups"
  };
  static const char *fn[3] = {
    "plus-plus.xvg",
    "min-min.xvg",
    "plus-min.xvg"
  };
  int        nset[3]={0,0,0};
  
  t_topology *top;
  int        ePBC;
  char       *buf;
  t_trxstatus *status;
  int        i,j,k,m,nnn,teller,ncg,n1,n2,n3,natoms;
  real       t,*time,qi,qj;
  t_charge   *cg;
  real       ***cgdist;
  int        **nWithin;
  
  double     t0,dt;
  char       label[234];
  t_pbc      pbc;
  rvec       *x;
  matrix     box;
  output_env_t oenv;
  
  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);
  cg=mk_charge(&top->atoms,&(top->cgs),&ncg);
  snew(cgdist,ncg);
  snew(nWithin,ncg);
  for(i=0; (i<ncg); i++) {
    snew(cgdist[i],ncg);
    snew(nWithin[i],ncg);
  }
  
  natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  
  teller=0;
  time=NULL;
  do {
    srenew(time,teller+1);
    time[teller]=t;
    
    set_pbc(&pbc,ePBC,box);
    
    for(i=0; (i<ncg); i++) {
      for(j=i+1; (j<ncg); j++) {
	srenew(cgdist[i][j],teller+1);
	cgdist[i][j][teller]=
	  calc_dist(&pbc,x,&(top->cgs),cg[i].cg,cg[j].cg);
	if (cgdist[i][j][teller] < truncate)
	  nWithin[i][j]=1;
      }
    }
          
    teller++;
  } while (read_next_x(oenv,status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);

  if (bSep) {
    snew(buf,256);
    for(i=0; (i<ncg); i++)
      for(j=i+1; (j<ncg); j++) {
	if (nWithin[i][j]) {
	  sprintf(buf,"sb-%s:%s.xvg",cg[i].label,cg[j].label);
	  fp=xvgropen(buf,buf,"Time (ps)","Distance (nm)",oenv);
	  for(k=0; (k<teller); k++) 
	    fprintf(fp,"%10g  %10g\n",time[k],cgdist[i][j][k]);
	  ffclose(fp);
	}
      }
    sfree(buf);
  }
  else {
  
    for(m=0; (m<3); m++)
      out[m]=xvgropen(fn[m],title[m],"Time (ps)","Distance (nm)",oenv);

    snew(buf,256);
    for(i=0; (i<ncg); i++) {
      qi=cg[i].q;
      for(j=i+1; (j<ncg); j++) {
	qj=cg[j].q;
	if (nWithin[i][j]) {
	  sprintf(buf,"%s:%s",cg[i].label,cg[j].label);
	  if (qi*qj < 0) 
	    nnn=2;
	  else if (qi+qj > 0) 
	    nnn=0;
	  else 
	    nnn=1;
	  
	  if (nset[nnn] == 0) 
	    xvgr_legend(out[nnn],1,(const char**)&buf,oenv);
	  else {
	    if (output_env_get_xvg_format(oenv) == exvgXMGR) {
	      fprintf(out[nnn],"@ legend string %d \"%s\"\n",nset[nnn],buf);
	    } else if (output_env_get_xvg_format(oenv) == exvgXMGRACE) {
	      fprintf(out[nnn],"@ s%d legend \"%s\"\n",nset[nnn],buf);
	    }
	  }
	  nset[nnn]++;
	  nWithin[i][j]=nnn+1;
	}
      }  
    }
    for(k=0; (k<teller); k++) {
      for(m=0; (m<3); m++)
	fprintf(out[m],"%10g",time[k]);
    
      for(i=0; (i<ncg); i++) {
	for(j=i+1; (j<ncg); j++) {
	  nnn=nWithin[i][j];
	  if (nnn >0) 
	    fprintf(out[nnn-1],"  %10g",cgdist[i][j][k]);
	}
      }
      for(m=0; (m<3); m++) 
	fprintf(out[m],"\n");
    }
    for(m=0; (m<3); m++) {
      ffclose(out[m]);
      if (nset[m] == 0)
	remove(fn[m]);
    }
  }
  thanx(stderr);
    
  return 0;
}
