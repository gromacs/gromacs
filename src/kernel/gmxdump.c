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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
static char *SRCID_gmxdump_c = "$Id$";

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
#include "assert.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "sheader.h"
#include "txtdump.h"

static void list_tpx(char *fn, bool bAltLayout)
{
  int         step,natoms,fp,indent;
  real        t,lambda;
  rvec        *x,*v,*f;
  matrix      box;
  t_inputrec  ir;
  t_tpxheader tpx;
  t_topology  top;
 
  read_tpxheader(fn,&tpx);
  snew(x,tpx.natoms);
  snew(v,tpx.natoms);
  snew(f,tpx.natoms);
  read_tpx(fn,&step,&t,&lambda,
	   tpx.bIr  ? &ir : NULL,
	   tpx.bBox ? box : NULL,
	   &natoms,
	   tpx.bX   ? x : NULL,
	   tpx.bV   ? v : NULL,
	   tpx.bF   ? f : NULL,
	   tpx.bTop ? &top: NULL);
  
  if (bAltLayout) {
    if (available(stdout,&tpx,fn))
      {
	indent=0;
	indent=pr_title(stdout,indent,fn);
	pr_header(stdout,indent,"header",&(tpx));
	pr_inputrec(stdout,indent,"ir",&(ir));
	pr_rvecs(stdout,indent,"box",box,DIM);
	pr_rvecs(stdout,indent,"x",x,natoms);
	pr_rvecs(stdout,indent,"v",v,natoms);
	pr_rvecs(stdout,indent,"f",f,natoms);
	pr_top(stdout,indent,"topology",&(top));
      }
  } else {
    printf("----------- begin of %s ------------\n",fn);
    fp = open_tpx(NULL,"w");
    fio_setdebug(fp,TRUE);
    fwrite_tpx(fp,step,t,lambda,
	       tpx.bIr  ? &ir : NULL,
	       tpx.bBox ? box : NULL,
	       natoms,
	       tpx.bX   ? x : NULL,
	       tpx.bV   ? v : NULL,
	       tpx.bF   ? f : NULL,
	       tpx.bTop ? &top : NULL);
    close_tpx(fp);
    printf("----------- end of %s ------------\n",fn);
  }
  
  sfree(x);
  sfree(v);
  sfree(f);
}

static void list_trn(char *fn)
{
  int         fpread,fpwrite,nframe;
  rvec        *x,*v,*f;
  matrix      box;
  t_trnheader trn;

  fpread  = open_trn(fn,"r"); 
  fpwrite = open_tpx(NULL,"w");
  fio_setdebug(fpwrite,TRUE);
  
  nframe = 0;
  while (fread_trnheader(fpread,&trn)) {
    snew(x,trn.natoms);
    snew(v,trn.natoms);
    snew(f,trn.natoms);
    fread_htrn(fpread,&trn,
	       trn.box_size ? box : NULL,
	       trn.x_size   ? x : NULL,
	       trn.v_size   ? v : NULL,
	       trn.f_size   ? f : NULL);
    printf("----------- begin of %s frame %d ------------\n",fn,nframe);
    fwrite_tpx(fpwrite,trn.step,trn.t,trn.lambda,NULL,
	       trn.box_size ? box : NULL,
	       trn.natoms,
	       trn.x_size   ? x : NULL,
	       trn.v_size   ? v : NULL,
	       trn.f_size   ? f : NULL,
	       NULL);
    printf("----------- end of %s frame %d ------------\n",fn,nframe);
    
    sfree(x);
    sfree(v);
    sfree(f);
    nframe++;
  }
  close_tpx(fpwrite);
  close_trn(fpread);
}

void list_xtc(char *fn)
{
  int    xd;
  rvec   *x;
  matrix box;
  int    natoms,step;
  real   prec,time;
  
  printf("gmxdump: %s\n",fn);
  
  xd = open_xtc(fn,"r");
  read_first_xtc(xd,&natoms,&step,&time,box,&x,&prec);
		 
  do {
    printf("natoms=%10d  step=%10d  time=%10g  prec=%10g\n",
	   natoms,step,time,prec);
    pr_rvecs(stdout,0,"box",box,DIM);
    pr_rvecs(stdout,0,"x",x,natoms);
  } while (read_next_xtc(xd,&natoms,&step,&time,box,x,&prec));
  
  close_xtc(xd);
}

void list_trx(char *fn)
{
  int ftp;
  
  ftp = fn2ftp(fn);
  if (ftp == efXTC)
    list_xtc(fn);
  else if ((ftp == efTRR) || (ftp == efTRJ))
    list_trn(fn);
  else
    fprintf(stderr,"File %s not supported. Try using more %s\n",
	    fn,fn);
}

void list_ene(char *fn,bool bEDR)
{
  int       in;
  bool      bCont;
  t_energy  *ener;
  t_drblock drblock;
  int       i,nre,nre2,step;
  real      t,rav,minthird;
  char      **enm=NULL;

  printf("gmxdump: %s\n",fn);
  in = open_enx(fn,"r");
  do_enxnms(in,&nre,&enm);
  
  printf("energy components:\n");
  for(i=0; (i<nre); i++) 
    printf("%5d  %s\n",i,enm[i]);
    
  drblock.ndr=0;
  minthird=-1.0/3.0;
  snew(ener,nre);
  do {
    bCont=do_enx(in,&t,&step,&nre2,ener,&drblock);
    assert(nre2==nre);
    
    if (bCont) {
      printf("\n%24s  %12.5e  %12s  %12d\n","time:",t,"step:",step);
      printf("%24s  %12s  %12s  %12s\n",
	     "Component","Energy","Av. Energy","Sum Energy");
      for(i=0; (i<nre); i++) 
	printf("%24s  %12.5e  %12.5e  %12.5e\n",
	       enm[i],ener[i].e,ener[i].eav,ener[i].esum);
      if (drblock.ndr)
	printf("Restraint %8s  %8s\n","r(t)","< r >");
      for(i=0; (i<drblock.ndr); i++) {
	rav=pow(drblock.rav[i],minthird);
	printf("%8d  %8.4f  %8.4f\n",i,drblock.rt[i],rav);
      }
    }
  } while (bCont);
  
  close_enx(in);
  
  sfree(ener);
  sfree(enm);
  if (drblock.ndr) {
    sfree(drblock.rav);
    sfree(drblock.rt);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gmxdump reads a run input file ([BB].tpa/,tpr/.tpb[bb]), a"
    "trajectory ([BB].trn[bb]) an energy ([BB].ene[bb] or [BB].edr[bb])",
    "or an xtc ([BB].xtc[bb]) file and",
    "prints that to standard output in a readable format.",
    "This program is essential for checking your run input file ",
    "in case of problems.[PAR]",
  };
  t_filenm fnm[] = {
    { efTPX, "-s", NULL, ffOPTRD },
    { efTRX, "-f", NULL, ffOPTRD },
    { efENX, "-e", NULL, ffOPTRD },
  };
#define NFILE asize(fnm)

  /* Command line options */
  static bool bAltLayout=FALSE;
  t_pargs pa[] = {
    { "-a", FALSE, etBOOL, &bAltLayout, "Alternative layout for run startup files" }
  };
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  if (ftp2bSet(efTPX,NFILE,fnm)) 
    list_tpx(ftp2fn(efTPX,NFILE,fnm), bAltLayout);
    
  if (ftp2bSet(efTRX,NFILE,fnm)) 
    list_trx(ftp2fn(efTRX,NFILE,fnm));
  
  if (ftp2bSet(efENX,NFILE,fnm))
    list_ene(ftp2fn(efENX,NFILE,fnm), FALSE);
    
  thanx(stdout);

  return 0;
}
