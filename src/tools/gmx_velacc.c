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
#include <stdio.h>
#include <math.h>

#include "confio.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "index.h"
#include "smalloc.h"
#include "statutil.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "strdb.h"
#include "xvgr.h"
#include "gmx_ana.h"
#include "gmx_fft.h"

static void index_atom2mol(int *n,atom_id *index,t_block *mols)
{
  int nat,i,nmol,mol,j;

  nat = *n;
  i = 0;
  nmol = 0;
  mol = 0;
  while (i < nat) {
    while (index[i] > mols->index[mol]) {
      mol++;
      if (mol >= mols->nr)
	gmx_fatal(FARGS,"Atom index out of range: %d",index[i]+1);
    }
    for(j=mols->index[mol]; j<mols->index[mol+1]; j++) {
      if (i >= nat || index[i] != j)
      	gmx_fatal(FARGS,"The index group does not consist of whole molecules");
      i++;
    }
    index[nmol++] = mol;
  }

  fprintf(stderr,"\nSplit group of %d atoms into %d molecules\n",nat,nmol);

  *n = nmol;
}

static void precalc(t_topology top,real normm[]){

  real mtot;
  int i,j,k,l;

  for(i=0;i<top.mols.nr;i++){
    k=top.mols.index[i];
    l=top.mols.index[i+1];
    mtot=0.0;

    for(j=k;j<l;j++)
      mtot+=top.atoms.atom[j].m;

    for(j=k;j<l;j++)
      normm[j]=top.atoms.atom[j].m/mtot;

  }

}

static void calc_spectrum(int n,real c[],real dt,const char *fn,
			  output_env_t oenv,gmx_bool bRecip)
{
  FILE *fp;
  gmx_fft_t fft;
  int  i,status;
  real *data;
  real nu,omega,recip_fac;

  snew(data,n*2);
  for(i=0; (i<n); i++)
    data[i] = c[i];

  if ((status = gmx_fft_init_1d_real(&fft,n,GMX_FFT_FLAG_NONE)) != 0)
    gmx_fatal(FARGS,"Invalid fft return status %d",status);

  if ((status = gmx_fft_1d_real(fft, GMX_FFT_REAL_TO_COMPLEX,data,data)) != 0)
    gmx_fatal(FARGS,"Invalid fft return status %d",status);

  fp = xvgropen(fn,"Vibrational Power Spectrum",
		bRecip ? "\\f{12}w\\f{4} (cm\\S-1\\N)" :
		"\\f{12}n\\f{4} (ps\\S-1\\N)",
		"a.u.",oenv);
  /* This is difficult.
   * The length of the ACF is dt (as passed to this routine).
   * We pass the vacf with N time steps from 0 to dt.
   * That means that after FFT we have lowest frequency = 1/dt
   * then 1/(2 dt) etc. (this is the X-axis of the data after FFT).
   * To convert to 1/cm we need to have to realize that
   * E = hbar w = h nu = h c/lambda. We want to have reciprokal cm
   * on the x-axis, that is 1/lambda, so we then have
   * 1/lambda = nu/c. Since nu has units of 1/ps and c has gromacs units
   * of nm/ps, we need to multiply by 1e7.
   * The timestep between saving the trajectory is
   * 1e7 is to convert nanometer to cm
   */
  recip_fac = bRecip ? (1e7/SPEED_OF_LIGHT) : 1.0;
  for(i=0; (i<n); i+=2) {
    nu = i/(2*dt);
    omega = nu*recip_fac;
    fprintf(fp,"%10g  %10g\n",omega,sqr(data[i])+sqr(data[i+1]));
  }
  xvgrclose(fp);
  gmx_fft_destroy(fft);
  sfree(data);
}

int gmx_velacc(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]g_velacc[tt] computes the velocity autocorrelation function.",
    "When the [TT]-m[tt] option is used, the momentum autocorrelation",
    "function is calculated.[PAR]",
    "With option [TT]-mol[tt] the velocity autocorrelation function of",
    "molecules is calculated. In this case the index group should consist",
    "of molecule numbers instead of atom numbers.[PAR]",
    "Be sure that your trajectory contains frames with velocity information",
    "(i.e. [TT]nstvout[tt] was set in your original [TT].mdp[tt] file),",
    "and that the time interval between data collection points is",
    "much shorter than the time scale of the autocorrelation."
  };
  
  static gmx_bool bM=FALSE,bMol=FALSE,bRecip=TRUE;
  t_pargs pa[] = {
    { "-m", FALSE, etBOOL, {&bM},
      "Calculate the momentum autocorrelation function" },
    { "-recip", FALSE, etBOOL, {&bRecip},
      "Use cm^-1 on X-axis instead of 1/ps for spectra." },
    { "-mol", FALSE, etBOOL, {&bMol},
      "Calculate the velocity acf of molecules" }
  };

  t_topology top;
  int        ePBC=-1;
  t_trxframe fr;
  matrix     box;
  gmx_bool       bTPS=FALSE,bTop=FALSE;
  int        gnx;
  atom_id    *index;
  char       *grpname;
  char       title[256];
  real       t0,t1,dt,m;
  t_trxstatus *status;
  int        teller,n_alloc,i,j,tel3,k,l;
  rvec       mv_mol;
  real       **c1;
  real	     *normm=NULL;
  output_env_t oenv;
  
#define NHISTO 360
    
  t_filenm  fnm[] = {
    { efTRN, "-f",    NULL,   ffREAD  },
    { efTPS, NULL,    NULL,   ffOPTRD }, 
    { efNDX, NULL,    NULL,   ffOPTRD },
    { efXVG, "-o",    "vac",  ffWRITE },
    { efXVG, "-os",   "spectrum", ffOPTWR }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL,&oenv);

  if (bMol || bM) {
    bTPS = ftp2bSet(efTPS,NFILE,fnm) || !ftp2bSet(efNDX,NFILE,fnm);
  }

  if (bTPS) {
    bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,NULL,NULL,box,
		       TRUE);
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  } else
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

  if (bMol) {
    if (!bTop)
      gmx_fatal(FARGS,"Need a topology to determine the molecules");
    snew(normm,top.atoms.nr);
    precalc(top,normm);
    index_atom2mol(&gnx,index,&top.mols);
  }
  
  /* Correlation stuff */
  snew(c1,gnx);
  for(i=0; (i<gnx); i++)
    c1[i]=NULL;
  
  read_first_frame(oenv,&status,ftp2fn(efTRN,NFILE,fnm),&fr,TRX_NEED_V);
  t0=fr.time;
      
  n_alloc=0;
  teller=0;
  do {
    if (teller >= n_alloc) {
      n_alloc+=100;
      for(i=0; i<gnx; i++)
	srenew(c1[i],DIM*n_alloc);
    }
    tel3=3*teller;
    if (bMol)
      for(i=0; i<gnx; i++) {
	clear_rvec(mv_mol);
	k=top.mols.index[index[i]];
	l=top.mols.index[index[i]+1];
	for(j=k; j<l; j++) {
   	  if (bM)
	    m = top.atoms.atom[j].m;
	  else
	    m = normm[j];
	  mv_mol[XX] += m*fr.v[j][XX];
	  mv_mol[YY] += m*fr.v[j][YY];
	  mv_mol[ZZ] += m*fr.v[j][ZZ];
	}
	c1[i][tel3+XX]=mv_mol[XX];
	c1[i][tel3+YY]=mv_mol[YY];
	c1[i][tel3+ZZ]=mv_mol[ZZ];
      }
     else
      for(i=0; i<gnx; i++) {
        if (bM)
	  m = top.atoms.atom[index[i]].m;
	else
	 m = 1;
	c1[i][tel3+XX]=m*fr.v[index[i]][XX];
	c1[i][tel3+YY]=m*fr.v[index[i]][YY];
	c1[i][tel3+ZZ]=m*fr.v[index[i]][ZZ];
      }

    t1=fr.time;

    teller ++;
  } while (read_next_frame(oenv,status,&fr));
  
  close_trj(status);
  dt = (t1-t0)/(teller-1);
  do_autocorr(opt2fn("-o",NFILE,fnm), oenv,
	      bM ? 
	      "Momentum Autocorrelation Function" :
	      "Velocity Autocorrelation Function",
	      teller,gnx,c1,dt,eacVector,TRUE);

  do_view(oenv,opt2fn("-o",NFILE,fnm),"-nxy");
  
  if (opt2bSet("-os",NFILE,fnm)) {
    calc_spectrum(teller/2,(real *) (c1[0]),(t1-t0)/2,opt2fn("-os",NFILE,fnm),
		  oenv,bRecip);
    do_view(oenv,opt2fn("-os",NFILE,fnm),"-nxy");
  }

  thanx(stderr);
  
  return 0;
}
