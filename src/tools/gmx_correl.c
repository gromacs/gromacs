/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */

/*
 * gmx_correl.c (c) 2010 Alexey Shvetsov 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <string.h>
#include <time.h>

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <direct.h>
#include <io.h>
#endif

#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "confio.h"
#include "trnio.h"
#include "mshift.h"
#include "xvgr.h"
#include "do_fit.h"
#include "rmpbc.h"
#include "txtdump.h"
#include "matio.h"
#include "physics.h"
#include "gmx_ana.h"

int gmx_correl(int argc,char *argv[]) {
  const char *desc[] = {
    "[TT]g_correl[tt] - simple tool to calculates and the correlation matrix.",
    "All structures are fitted to the structure in the structure file.",
    "When this is not a run input file periodicity will not be taken into",
    "account.",
    "[PAR]",
    "When the same atoms are used for the fit and the correlation analysis,",
    "the reference structure for the fit is written first with t=-1.",
    "[PAR]",
    "Option [TT]-ascii[tt] writes the whole correlation matrix to",
    "an ASCII file.", 
    "[PAR]",
    "Option [TT]-xpm[tt] writes the whole correlation matrix to an xpm file.",
    "[PAR]"
  };
  gmx_bool bFit=TRUE,bRef=FALSE,bPBC=TRUE;
  t_pargs pa[] = {
    { "-fit",  FALSE, etBOOL, {&bFit},
      "Fit to a reference structure"},
    { "-ref",  FALSE, etBOOL, {&bRef},
      "Use the deviation from the conformation in the structure file instead of from the average" },
    { "-pbc",  FALSE,  etBOOL, {&bPBC},
      "Apply corrections for periodic boundary conditions" }
  };
  FILE          *out;
  t_trxstatus   *status;
  t_topology    top;
  int           ePBC;
  t_atoms       *atoms;
  rvec          *xread,*xref;
  rvec          *xa,*xava;
  rvec          *xb,*xavb;
  matrix        box,zerobox;
  real          inv_nframes;
  real          t,tstart,tend;
  real          *w_rls=NULL;
  real          min,max,*axisa,*axisb;
  real          **mat;
  real          *r2a,*r2b;
  int           natoms,nat,nframes0,nframes,nlevels;
  int           natomsa,natomsb;
  gmx_large_int_t i,j;
  const char *fitfile,*trxfile,*ndxfile;
  const char *asciifile,*xpmfile,*averfilea,*averfileb;
  char       str[STRLEN],*fitname;
  char       *ananamea,*ananameb;
  int        d,nfit;
  atom_id    *ifit;
  atom_id    *indexa,*indexb;
  t_rgb      rlo,rmi,rhi;
  output_env_t oenv;
  gmx_rmpbc_t  gpbc=NULL;

  t_filenm fnm[] = { 
    { efTRX, "-f",  NULL, ffREAD }, 
    { efTPS, NULL,  NULL, ffREAD },
    { efNDX, NULL,  NULL, ffOPTRD },
    { efSTO, "-ava", "averagea.pdb", ffOPTWR },
    { efSTO, "-avb", "averageb.pdb", ffOPTWR },
    { efDAT, "-ascii","correlation", ffWRITE },
    { efXPM, "-xpm","correlation", ffOPTWR }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv); 

  clear_mat(zerobox);

  fitfile    = ftp2fn(efTPS,NFILE,fnm);
  trxfile    = ftp2fn(efTRX,NFILE,fnm);
  ndxfile    = ftp2fn_null(efNDX,NFILE,fnm);
  averfilea  = ftp2fn(efSTO,NFILE,fnm);
  averfileb  = ftp2fn(efSTO,NFILE,fnm);
  asciifile  = opt2fn_null("-ascii",NFILE,fnm);
  xpmfile    = opt2fn_null("-xpm",NFILE,fnm);

  read_tps_conf(fitfile,str,&top,&ePBC,&xref,NULL,box,TRUE);
  atoms=&top.atoms;

  if (bFit) {
    printf("\nChoose a group for the least squares fit\n"); 
    get_index(atoms,ndxfile,1,&nfit,&ifit,&fitname);
    if (nfit < 3) 
      gmx_fatal(FARGS,"Need >= 3 points to fit!\n");
  } else
    nfit=0;
  printf("\nChoose two groups for the correlation analysis\n"); 
  get_index(atoms,ndxfile,1,&natomsa,&indexa,&ananamea);
  get_index(atoms,ndxfile,1,&natomsb,&indexb,&ananameb);

  if (bFit) {
    snew(w_rls,atoms->nr);
    for(i=0; (i<nfit); i++)
	w_rls[ifit[i]]=1.0;
    }

  /* Prepare reference frame */
  if (bPBC)
   gpbc = gmx_rmpbc_init(&top.idef,ePBC,atoms->nr,box);
   gmx_rmpbc(gpbc,atoms->nr,box,xref);
  if (bFit)
   reset_x(nfit,ifit,atoms->nr,NULL,xref,w_rls);

  snew(xa,natomsa);
  snew(xb,natomsb);
  snew(xava,natomsa);
  snew(xavb,natomsb);

  fprintf(stderr,"Calculating the average structure ...\n");
  nframes0 = 0;
  nat=read_first_x(oenv,&status,trxfile,&t,&xread,box);
  if (nat != atoms->nr)
    fprintf(stderr,"\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n",natoms,nat);
  do {
    nframes0++;
    /* calculate x: a fitted struture of the selected atoms */
    if (bPBC)
      gmx_rmpbc(gpbc,atoms->nr,box,xread);
    if (bFit) {
      reset_x(nfit,ifit,nat,NULL,xread,w_rls);
      do_fit(nat,w_rls,xref,xread);
    }
    for (i=0; i<natomsa; i++)
      rvec_inc(xava[i],xread[indexa[i]]);
    for (i=0; i<natomsb; i++)
      rvec_inc(xavb[i],xread[indexb[i]]);
  } while (read_next_x(oenv,status,&t,nat,xread,box));
  close_trj(status);
  
  inv_nframes = 1.0/nframes0;
  for(i=0; i<natomsa; i++)
    for(d=0; d<DIM; d++) {
      xava[i][d] *= inv_nframes;
      xread[indexa[i]][d] = xava[i][d];
    }
  write_sto_conf_indexed(opt2fn("-ava",NFILE,fnm),"Average structure",
  		         atoms,xread,NULL,epbcNONE,zerobox,natomsa,indexa);
  for(i=0; i<natomsb; i++)
   for(d=0; d<DIM; d++) {
      xavb[i][d] *= inv_nframes;
      xread[indexb[i]][d] = xavb[i][d];
   }
  write_sto_conf_indexed(opt2fn("-avb",NFILE,fnm),"Average structure",
			 atoms,xread,NULL,epbcNONE,zerobox,natomsb,indexb);
  sfree(xread);

  fprintf(stderr,"Constructing correlation matrix (%dx%d) ...\n",(int)natomsa,(int)natomsb);
  nframes=0;
  nat=read_first_x(oenv,&status,trxfile,&t,&xread,box);
  tstart = t;
  snew(r2a,natomsa);
  snew(r2b,natomsb);
  /* Allocate 2d matrix */
  snew(mat,natomsa);
  for (i=0;i<natomsb;i++)
  	snew(mat[i],natomsb);
  do {
    nframes++;
    tend = t;
    /* calculate x: a (fitted) structure of the selected atoms */
    if (bPBC)
      gmx_rmpbc(gpbc,atoms->nr,box,xread);
    if (bFit) {
      reset_x(nfit,ifit,nat,NULL,xread,w_rls);
      do_fit(nat,w_rls,xref,xread);
    }
    if (bRef) {
      for (i=0; i<natomsa; i++)
	rvec_sub(xread[indexa[i]],xref[indexa[i]],xa[i]);
      for (i=0; i<natomsb;i++)
        rvec_sub(xread[indexb[i]],xref[indexb[i]],xb[i]);
      }
    else {
      for (i=0; i<natomsa; i++)
	rvec_sub(xread[indexa[i]],xava[i],xa[i]);
      for (i=0; i<natomsb; i++)
        rvec_sub(xread[indexb[i]],xavb[i],xb[i]);
      }
    /* find actual rmsf */
    for (i=0; i<natomsa; i++)
        r2a[i]+=norm2(xa[i]);
    for (i=0; i<natomsb; i++)
        r2b[i]+=norm2(xb[i]);
    /* TODO: Optimize calculation of correlation matrix */
    for (i=0; i<natomsa; i++) {
    	for (j=0; j<natomsb; j++) {
		mat[i][j] += iprod(xa[i],xb[j]);
	}
    }
  } while (read_next_x(oenv,status,&t,nat,xread,box) && 
	   (bRef || nframes < nframes0));
  close_trj(status);
  gmx_rmpbc_done(gpbc);
  fprintf(stderr,"Read %d frames\n",nframes);
  /* Correct mat so it will be real correlation one */
  for (i=0;i<natomsa;i++) {
  	for (j=0;j<natomsb;j++) {
		mat[i][j] = mat[i][j]/sqrt(r2a[i]*r2b[j]);
	}
  }
  if (asciifile) {
    out = ffopen(asciifile,"w");
    for (i=0; i<natomsa; i++) {
      for (j=0; j<natomsb; j++) {
	fprintf(out,"%g ",mat[i][j]);
	}
      fprintf(out,"\n");
    }
    ffclose(out);
  }
  
  if (xpmfile) {
    min = 0;
    max = 0;
    for (i=0; i<natomsa; i++) {
      for (j=0; j<natomsb; j++) {
	if (mat[i][j] < min)
	  min = mat[i][j];
	if (mat[i][j] > max)
	  max = mat[i][j];
      }
    }
    snew(axisa,natomsa);
    snew(axisb,natomsb);
    for(i=0;i<natomsa;i++)
      axisa[i] = i+1;
    for(i=0;i<natomsb;i++)
      axisb[i] = i+1;
    rlo.r = 0; rlo.g = 0; rlo.b = 1;
    rmi.r = 1; rmi.g = 1; rmi.b = 1;
    rhi.r = 1; rhi.g = 0; rhi.b = 0;
    out = ffopen(xpmfile,"w");
    nlevels = 80;
    write_xpm3(out,0,"Correlation", "",
	       "Residue Number","Residue Number",natomsa,natomsb,axisa,axisb,
	       mat,min,0.0,max,rlo,rmi,rhi,&nlevels);
    ffclose(out);
    sfree(axisa);
    sfree(axisb);
  }
  for (i=0;i<natomsb;i++)
  	sfree(mat[i]);
  sfree(mat);
  sfree(r2a);
  sfree(r2b);

  thanx(stderr);

  return 0;
}
