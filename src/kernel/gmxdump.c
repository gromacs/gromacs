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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "macros.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "fatal.h"
#include "xtcio.h"
#include "enxio.h"
#include "smalloc.h"
#include "names.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "txtdump.h"

static void list_tpx(char *fn, bool bShowNumbers)
{
  int         step,fp,indent,i,j,**gcount,atot;
  real        t;
  t_state     state;
  rvec        *f=NULL;
  t_inputrec  ir;
  t_tpxheader tpx;
  t_topology  top;
  read_tpxheader(fn,&tpx,TRUE,NULL,NULL);

  read_tpx_state(fn,&step,&t,
		 tpx.bIr  ? &ir : NULL,
		 &state,tpx.bF ? f : NULL,
		 tpx.bTop ? &top: NULL);
  
  if (available(stdout,&tpx,fn)) {
    indent=0;
    indent=pr_title(stdout,indent,fn);
    pr_header(stdout,indent,"header",&(tpx));
    pr_inputrec(stdout,indent,"ir",&(ir));
    pr_rvecs(stdout,indent,"box",tpx.bBox ? state.box : NULL,DIM);
    pr_rvecs(stdout,indent,"boxv",tpx.bBox ? state.boxv : NULL,DIM);
    pr_rvecs(stdout,indent,"pcoupl_mu",tpx.bBox ? state.pcoupl_mu : NULL,DIM);
    pr_reals(stdout,indent,"nosehoover_xi",state.nosehoover_xi,state.ngtc);
    pr_reals(stdout,indent,"tcoupl_lambda",state.tcoupl_lambda,state.ngtc);
    pr_rvecs(stdout,indent,"x",tpx.bX ? state.x : NULL,state.natoms);
    pr_rvecs(stdout,indent,"v",tpx.bV ? state.v : NULL,state.natoms);
    pr_rvecs(stdout,indent,"f",f,state.natoms);
    pr_top(stdout,indent,"topology",&(top),bShowNumbers);
  }

  snew(gcount,egcNR);
  for(i=0; (i<egcNR); i++) 
    snew(gcount[i],top.atoms.grps[i].nr);
  
  for(i=0; (i<top.atoms.nr); i++) {
    for(j=0; (j<egcNR); j++) 
      gcount[j][top.atoms.atom[i].grpnr[j]]++;
  }
  printf("Group statistics\n");
  for(i=0; (i<egcNR); i++) {
    atot=0;
    printf("%-12s: ",gtypes[i]);
    for(j=0; (j<top.atoms.grps[i].nr); j++) {
      printf("  %5d",gcount[i][j]);
      atot+=gcount[i][j];
    }
    printf("  (total %d atoms)\n",atot);
    sfree(gcount[i]);
  }
  sfree(gcount);

  done_state(&state);
  sfree(f);
}

static void list_trn(char *fn)
{
  int         fpread,fpwrite,nframe,indent;
  char        buf[256];
  rvec        *x,*v,*f;
  matrix      box;
  t_trnheader trn;
  bool        bOK;

  fpread  = open_trn(fn,"r"); 
  fpwrite = open_tpx(NULL,"w");
  fio_setdebug(fpwrite,TRUE);
  
  nframe = 0;
  while (fread_trnheader(fpread,&trn,&bOK)) {
    snew(x,trn.natoms);
    snew(v,trn.natoms);
    snew(f,trn.natoms);
    if (fread_htrn(fpread,&trn,
		   trn.box_size ? box : NULL,
		   trn.x_size   ? x : NULL,
		   trn.v_size   ? v : NULL,
		   trn.f_size   ? f : NULL)) {
      sprintf(buf,"%s frame %d",fn,nframe);
      indent=0;
      indent=pr_title(stdout,indent,buf);
      pr_indent(stdout,indent);
      fprintf(stdout,"natoms=%10d  step=%10d  time=%10g  lambda=%10g\n",
	      trn.natoms,trn.step,trn.t,trn.lambda);
      if (trn.box_size)
	pr_rvecs(stdout,indent,"box",box,DIM);
      if (trn.x_size)
	pr_rvecs(stdout,indent,"x",x,trn.natoms);
      if (trn.v_size)
	pr_rvecs(stdout,indent,"v",v,trn.natoms);
      if (trn.f_size)
	pr_rvecs(stdout,indent,"f",f,trn.natoms);
    } 
    else
      fprintf(stderr,"\nWARNING: Incomplete frame: nr %d, t=%g\n",
	      nframe,trn.t);
    
    sfree(x);
    sfree(v);
    sfree(f);
    nframe++;
  }
  if (!bOK)
    fprintf(stderr,"\nWARNING: Incomplete frame header: nr %d, t=%g\n",
	    nframe,trn.t);
  close_tpx(fpwrite);
  close_trn(fpread);
}

void list_xtc(char *fn, bool bXVG)
{
  int    xd,indent;
  char   buf[256];
  rvec   *x;
  matrix box;
  int    nframe,natoms,step;
  real   prec,time;
  bool   bOK;
  
  xd = open_xtc(fn,"r");
  read_first_xtc(xd,&natoms,&step,&time,box,&x,&prec,&bOK);
		
  nframe=0;
  do {
    if (bXVG) {
      int i,d;
      
      fprintf(stdout,"%g",time);
      for(i=0; i<natoms; i++)
	for(d=0; d<DIM; d++)
	  fprintf(stdout," %g",x[i][d]);
      fprintf(stdout,"\n");
    } else {
      sprintf(buf,"%s frame %d",fn,nframe);
      indent=0;
      indent=pr_title(stdout,indent,buf);
      pr_indent(stdout,indent);
      fprintf(stdout,"natoms=%10d  step=%10d  time=%10g  prec=%10g\n",
	    natoms,step,time,prec);
      pr_rvecs(stdout,indent,"box",box,DIM);
      pr_rvecs(stdout,indent,"x",x,natoms);
    }
    nframe++;
  } while (read_next_xtc(xd,natoms,&step,&time,box,x,&prec,&bOK));
  if (!bOK)
    fprintf(stderr,"\nWARNING: Incomplete frame at time %g\n",time);
  close_xtc(xd);
}

void list_trx(char *fn,bool bXVG)
{
  int ftp;
  
  ftp = fn2ftp(fn);
  if (ftp == efXTC)
    list_xtc(fn,bXVG);
  else if ((ftp == efTRR) || (ftp == efTRJ))
    list_trn(fn);
  else
    fprintf(stderr,"File %s not supported. Try using more %s\n",
	    fn,fn);
}

void list_ene(char *fn)
{
  int        in,ndr;
  bool       bCont;
  t_enxframe *fr;
  int        i,nre,b;
  real       rav,minthird;
  char       **enm=NULL;

  printf("gmxdump: %s\n",fn);
  in = open_enx(fn,"r");
  do_enxnms(in,&nre,&enm);
  
  printf("energy components:\n");
  for(i=0; (i<nre); i++) 
    printf("%5d  %s\n",i,enm[i]);
    
  minthird=-1.0/3.0;
  snew(fr,1);
  do {
    bCont=do_enx(in,fr);
    if (fr->nre != nre)
      gmx_incons("Number of energies in list_ene");
    
    if (bCont) {
      printf("\n%24s  %12.5e  %12s  %12d\n","time:",
	     fr->t,"step:",fr->step);
      printf("%24s  %12s  %12s  %12s\n",
	     "Component","Energy","Av. Energy","Sum Energy");
      for(i=0; (i<nre); i++) 
	printf("%24s  %12.5e  %12.5e  %12.5e\n",
	       enm[i],fr->ener[i].e,fr->ener[i].eav,fr->ener[i].esum);
      if (fr->ndisre > 0) {
	printf("Distance restraint %8s  %8s\n","r(t)","< r >");
	for(i=0; i<fr->ndisre; i++) {
	  rav=pow(fr->rav[i],minthird);
	  printf("%17d  %8.4f  %8.4f\n",i,fr->rt[i],rav);
	}
      }
      for(b=0; b<fr->nblock; b++)
	if (fr->nr[b] > 0) {
	  printf("Block data %2d (%4d elm.) %8s\n",b,fr->nr[b],"value");
	  for(i=0; i<fr->nr[b]; i++)
	    printf("%24d  %8.4f\n",i,fr->block[b][i]);
	}
    }
  } while (bCont);
  
  close_enx(in);

  free_enxframe(fr);
  sfree(fr);
  sfree(enm);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gmxdump reads a run input file ([TT].tpa[tt]/[TT].tpr[tt]/[TT].tpb[tt]),",
    "a trajectory ([TT].trj[tt]/[TT].trr[tt]/[TT].xtc[tt]) or an energy",
    "file ([TT].ene[tt]/[TT].edr[tt]) and prints that to standard",
    "output in a readable format. This program is essential for",
    "checking your run input file in case of problems.[PAR]"
  };
  t_filenm fnm[] = {
    { efTPX, "-s", NULL, ffOPTRD },
    { efTRX, "-f", NULL, ffOPTRD },
    { efENX, "-e", NULL, ffOPTRD },
  };
#define NFILE asize(fnm)

  /* Command line options */
  static bool bXVG=FALSE;
  static bool bShowNumbers=TRUE;
  t_pargs pa[] = {
    { "-xvg", FALSE, etBOOL, {&bXVG}, "HIDDENXVG layout for xtc" },
    { "-nr",FALSE, etBOOL, {&bShowNumbers},"Show index numbers in output (leaving them out makes comparison easier, but creates a useless topology)" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);


  if (ftp2bSet(efTPX,NFILE,fnm)) 
    list_tpx(ftp2fn(efTPX,NFILE,fnm),bShowNumbers);
    
  if (ftp2bSet(efTRX,NFILE,fnm)) 
    list_trx(ftp2fn(efTRX,NFILE,fnm),bXVG);
  
  if (ftp2bSet(efENX,NFILE,fnm))
    list_ene(ftp2fn(efENX,NFILE,fnm));
    
  thanx(stderr);

  return 0;
}
