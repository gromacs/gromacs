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
static char *SRCID_gmxcheck_c = "$Id$";

#include <stdio.h>
#include <string.h>
#include "main.h"
#include "macros.h"
#include "math.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "fatal.h"
#include "sheader.h"
#include "statusio.h"
#include "xtcio.h"
#include "tpbcmp.h"

void chk_trj(int ftp,char *fn)
{
  t_statheader sh,count;
  int          idum,j,natoms,step;
  real         rdum,t,t0,old_t1,old_t2,prec;
  bool         bShowTimestep=TRUE;
  rvec         *x;
  matrix       box;
  size_t       fpos;
  XDR          xd;
  FILE         *status;
  
#define BOGUSTIME -1e10

  natoms = 0;  
  t      = 0;
  t0     = BOGUSTIME;
  step   = 0;
  count.e_size=0;
  count.box_size=0;
  count.vir_size=0;
  count.pres_size=0;
  count.x_size=0;
  count.v_size=0;
  count.f_size=0;

  printf("Checking file %s\n",fn);
  switch (ftp) {
  case efTRJ:
    status = ffopen(fn,"r");
    j=0;
    t=-1;
    old_t2=-2.0;
    old_t1=-1.0;
    while (!eof(status)) {
      fpos=ftell(status);
      rd_header(status,&sh);
      if (j>=2) {
	if ( fabs((sh.t-old_t1)-(old_t1-old_t2)) > 
	     0.1*(fabs(sh.t-old_t1)+fabs(old_t1-old_t2)) ) {
	  bShowTimestep=FALSE;
	  fprintf(stderr,"\nTimesteps at t=%g don't match (%g, %g)\n",
		  old_t1,old_t1-old_t2,sh.t-old_t1);
	}
      }
      old_t2=old_t1;
      old_t1=sh.t;
      if (t0 == BOGUSTIME) t0=sh.t;
      fprintf(stderr,"\rframe: %6d, t: %10.3f bytes: %10d",j,sh.t,fpos);
      if (j == 0)
	fprintf(stderr,"\n");
      rd_hstatus(status,&sh,&idum,&t,&rdum,NULL,
		 NULL,NULL,NULL,&natoms,NULL,NULL,NULL,&idum,NULL,NULL);
      j++;
#define INC(s,n,item) if (s.item  != 0) n.item++
      INC(sh,count,e_size);
      INC(sh,count,box_size);
      INC(sh,count,vir_size);
      INC(sh,count,pres_size);
      INC(sh,count,x_size);
      INC(sh,count,v_size);
      INC(sh,count,f_size);
    }
    fclose(status);
    t=sh.t;
    break;
  case efXTC:
    if (read_first_xtc(&xd,fn,&natoms,&step,&t,box,&x,&prec)) {
      fprintf(stderr,"\nXTC precision %g\n\n",prec);
      j=0;
      old_t2=-2.0;
      old_t1=-1.0;
      do {
	if (j>=2) {
	  if ( fabs((t-old_t1)-(old_t1-old_t2)) > 
	       0.1*(fabs(t-old_t1)+fabs(old_t1-old_t2)) ) {
	    bShowTimestep=FALSE;
	    fprintf(stderr,"\nTimesteps at t=%g don't match (%g, %g)\n",
		    old_t1,old_t1-old_t2,t-old_t1);
	  }
	}
	old_t2=old_t1;
	old_t1=t;
	if (t0 == BOGUSTIME) t0=t;
	fprintf(stderr,"\rframe: %6d, t: %10.3f",j,t);
	if (j == 0)
	  fprintf(stderr,"\n");
	count.x_size++;
	count.box_size++;
	j++;
      } while (read_next_xtc(&xd,&natoms,&step,&t,box,x,&prec));
      close_xtc(&xd);
    }
    else
      fprintf(stderr,"Empty file %s\n",fn);
    break;
  default:
    fprintf(stderr,"Sorry %s not supported yet\n",fn);
  }

  fprintf(stderr,"\n\n");
  fprintf(stderr,"\n# Atoms     %d\n",natoms);
  /*
  if (bShowTimestep) {
    fprintf(stderr,"\nItem        Timestep (ps)\n");
    if ((bShowTimestep) && (count.x_size > 1))
      fprintf(stderr,"Coords      %g\n",(t-t0)/(count.x_size-1));
    if (count.v_size > 1)
      fprintf(stderr,"Velocities  %g\n",(t-t0)/(count.v_size-1));
    if (count.f_size > 1)
      fprintf(stderr,"Forces      %g\n",(t-t0)/(count.f_size-1));
  }
  */
  fprintf(stderr,"\nItem        #frames");
  if (bShowTimestep)
    fprintf(stderr," Timestep (ps)");
  fprintf(stderr,"\n");
#define PRINTITEM(label,item) fprintf(stderr,"%-10s  %6d",label,count.item); if ((bShowTimestep) && (count.item > 1)) fprintf(stderr,"    %g\n",(t-t0)/(count.item-1)); else fprintf(stderr,"\n")
  PRINTITEM ( "Energies",   e_size );
  PRINTITEM ( "Box",        box_size );
  PRINTITEM ( "Virial",     vir_size );
  PRINTITEM ( "Pressure",   pres_size );
  PRINTITEM ( "Coords",     x_size );
  PRINTITEM ( "Velocities", v_size );
  PRINTITEM ( "Forces",     f_size );
  /*
    fprintf(stderr,"Energies    %d\n",count.e_size);
    fprintf(stderr,"Box         %d\n",count.box_size);
    fprintf(stderr,"Virial      %d\n",count.vir_size);
    fprintf(stderr,"Pressure    %d\n",count.pres_size);
    fprintf(stderr,"Coords      %d\n",count.x_size);
    fprintf(stderr,"Velocities  %d\n",count.v_size);
    fprintf(stderr,"Forces      %d\n",count.f_size);
    */
}  

void chk_ndx(char *fn)
{
  FILE *in;
  char buf[256];
  int i,j,ngrp,nelem,nitem,nat,ng_cnt,na_cnt;
  int nLab;

  fprintf(stderr,"Checking index file %s\n",fn);
  in=ffopen(fn,"r");
  
  if (fscanf(in,"%d%d",&ngrp,&nat) != 2) {
    fprintf(stderr,"Couldn't read NGRP or NAT in file %s\n",fn);
    exit(1);
  }
  nLab=2;
  ng_cnt=0,na_cnt=0;
  fprintf(stderr,"There should be %d groups containing %d atoms in total\n",
	  ngrp,nat);
  while (fscanf(in,"%s",buf) == 1) {
    if (nLab == 2) {
      fprintf(stderr,"Group %4d: %16s",ng_cnt,buf);
      ng_cnt++;
      nLab--;
    }
    else if (nLab == 1) {
      fprintf(stderr,"  Nelem: %16s\n",buf);
      if (sscanf(buf,"%d",&nelem) != 1) {
	fprintf(stderr,"\nERROR in indexfile %s\n",fn);
	fprintf(stderr,"Couldn't find proper NATOMS in buf '%s'\n",buf);
	fprintf(stderr,"May be you have entered too many atoms in the previous group\n\n");
	ffclose(in);
	return;
      }
      nLab--;
    }
    else {
      if (sscanf(buf,"%d",&nitem) != 1) {
	fprintf(stderr,"\nERROR in indexfile %s\n",fn);
	fprintf(stderr,"Couldn't find proper ATOM in buf '%s'\n",buf);
	fprintf(stderr,"You should have entered %d more atom(s) for the previous group\n\n",nelem);
	ffclose(in);
	return;
      }
      nelem--;
      na_cnt++;
      if (nelem == 0)
	nLab=2;
    }
  }
  fprintf(stderr,"\nFound %6d groups, should be %6d\n",ng_cnt,ngrp);
  fprintf(stderr,"Found %6d atoms,  should be %6d\n",na_cnt,nat);
  
  ffclose(in);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gmxcheck reads a binary trajectory ([TT].trj[tt]), ",
    "a xtc ([TT].xtc[tt]) or an index ([TT].ndx[tt]) file and",
    "prints out useful information about these files",
    "and their contents of course.[PAR]",
    "The program will compare binary topology ([TT].tpb[tt]) files",
    "when both [TT]-s1[tt] and [TT]-s2[tt] are supplied."
  };
  t_filenm fnm[] = {
    { efTRJ, "-f", NULL, ffOPTRD },
    { efXTC, "-x", NULL, ffOPTRD },
    { efNDX, "-n", NULL, ffOPTRD },
    { efTPB, "-s1", "top1", ffOPTRD },
    { efTPB, "-s2", "top2", ffOPTRD }

  };
#define NFILE asize(fnm)
  char *fn1=NULL,*fn2=NULL;

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,
		    asize(desc),desc,0,NULL);
  
  if (ftp2bSet(efTRJ,NFILE,fnm)) {
    chk_trj(efTRJ,ftp2fn(efTRJ,NFILE,fnm));
  }
  if (ftp2bSet(efXTC,NFILE,fnm)) {
    chk_trj(efXTC,ftp2fn(efXTC,NFILE,fnm));
  }
  if (ftp2bSet(efNDX,NFILE,fnm)) {
    chk_ndx(ftp2fn(efNDX,NFILE,fnm));
  }
  if (opt2bSet("-s1",NFILE,fnm))
    fn1=opt2fn("-s1",NFILE,fnm);
  if (opt2bSet("-s2",NFILE,fnm))
    fn2=opt2fn("-s2",NFILE,fnm);
  if (fn1 && fn2)
    comp_tpb(fn1,fn2);
  else if (fn1 || fn2)
    fprintf(stderr,"Please give me TWO .tpb files!\n");
  
  thanx(stderr);
  
  return 0;
}
