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
#include "tconf.h"
#include "xtcio.h"
#include "enerio.h"
#include "assert.h"
#include "smalloc.h"

static void list_file(char *fn)
{
  FILE     *fp;
  char     *version;
  t_config cf;

  printf("gmxdump: %s\n",fn);
  fp=ffopen(fn,"r");
  version=init_configuration(fp,&cf);
  fprintf(stderr,"version=%s\n",version);
  do {
    version=rd_configuration(fp,&cf); 
    pr_configuration(stdout,0,fn,&cf.data);
  } while (!eof(fp)); 
  ffclose(fp);
}

void list_xtc(char *fn)
{
  XDR    xd;
  rvec   *x;
  matrix box;
  int    natoms,step;
  real   prec,time;
  
  printf("gmxdump: %s\n",fn);
  
  read_first_xtc(&xd,fn,&natoms,&step,&time,box,&x,&prec);
		 
  do {
    printf("natoms=%10d  step=%10d  time=%10g  prec=%10g\n",
	   natoms,step,time,prec);
    pr_rvecs(stdout,0,"box",box,DIM);
    pr_rvecs(stdout,0,"x",x,natoms);
  } while (read_next_xtc(&xd,&natoms,&step,&time,box,x,&prec));
  
  close_xtc(&xd);
}

void list_ene(char *fn,bool bEDR)
{
  FILE      *in;
  XDR       xdrs;
  bool      bCont;
  t_energy  *ener;
  t_drblock drblock;
  int       i,nre,nre2,step;
  real      t,rav,minthird;
  char      **enm=NULL;

  printf("gmxdump: %s\n",fn);
  if (bEDR) {
    xdropen(&xdrs,fn,"r");
    edr_nms(&xdrs,&nre,&enm);
  }
  else {
    in=ffopen(fn,"r");
    rd_ener_nms(in,&nre,&enm);
  }
  printf("energy components:\n");
  for(i=0; (i<nre); i++) 
    printf("%5d  %s\n",i,enm[i]);
    
  drblock.ndr=0;
  minthird=-1.0/3.0;
  snew(ener,nre);
  do {
    if (bEDR) {
      bCont=edr_io(&xdrs,&t,&step,&nre2,ener,&drblock);
      assert(nre2==nre);
    }
    else
      bCont=rd_ener(in,&t,&step,ener,&drblock);
      
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
  
  if (bEDR)
    xdrclose(&xdrs);
  else
    ffclose(in);
  
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
    "rstat reads a binary topology ([BB].tpb[bb]), a"
    "trajectory ([BB].trj[bb]) an energy ([BB].ene[bb] or [BB].edr[bb])",
    " or a xtc ([BB].xtc[bb]) file and",
    "prints that to standard output in a readable format.",
    "This program is essential for checking your topology",
    "in case of problems.[PAR]",
  };
  t_filenm fnm[] = {
    { efTPB, "-s", NULL, ffOPTRD },
    { efXTC, "-x", NULL, ffOPTRD },
    { efTRJ, "-f", NULL, ffOPTRD },
    { efENE, "-e", NULL, ffOPTRD },
    { efEDR, "-d", NULL, ffOPTRD }
  };
#define NFILE asize(fnm)
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,
		    asize(desc),desc,0,NULL);
  
  if (ftp2bSet(efTPB,NFILE,fnm)) 
    list_file(ftp2fn(efTPB,NFILE,fnm));
    
  if (ftp2bSet(efTRJ,NFILE,fnm)) 
    list_file(ftp2fn(efTRJ,NFILE,fnm));
  
  if (ftp2bSet(efXTC,NFILE,fnm)) 
    list_xtc(ftp2fn(efXTC,NFILE,fnm));

  if (ftp2bSet(efENE,NFILE,fnm))
    list_ene(ftp2fn(efENE,NFILE,fnm),FALSE);
    
  if (ftp2bSet(efEDR,NFILE,fnm))
    list_ene(ftp2fn(efEDR,NFILE,fnm),TRUE);
    
  thanx(stdout);

  return 0;
}
