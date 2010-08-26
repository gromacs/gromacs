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
#include "pbc.h"
#include "gmx_ana.h"


#define NK  24
#define NPK 4

#define NKC  6
#define NKC0 4
int kset_c[NKC+1] = { 0, 3, 9, 13, 16, 19, NK };

rvec v0[NK]={{1,0,0},{0,1,0},{0,0,1}, {1, 1,0},{1,-1,0},{1,0,1}, {1,0,-1},{0,1, 1},{0,1,-1}, {1, 1, 1},{1, 1,-1},{1,-1, 1},{-1,1, 1}, {2,0,0},{0,2,0},{0,0,2}, {3,0,0},{0,3,0},{0,0,3}, {4,0,0},{0,4,0},{0,0,4}};
rvec v1[NK]={{0,1,0},{0,0,1},{1,0,0}, {0, 0,1},{0, 0,1},{0,1,0}, {0,1, 0},{1,0, 0},{1,0, 0}, {1,-1, 0},{1,-1, 0},{1, 0,-1},{ 0,1,-1}, {0,1,0},{0,0,1},{1,0,0}, {0,1,0},{0,0,1},{1,0,0}, {0,1,0},{0,0,1},{1,0,0}};
rvec v2[NK]={{0,0,1},{1,0,0},{0,1,0}, {1,-1,0},{1, 1,0},{1,0,-1},{1,0, 1},{0,1,-1},{0,1, 1}, {1, 1,-2},{1, 1, 2},{1, 2, 1},{ 2,1, 1}, {0,0,1},{1,0,0},{0,1,0}, {0,0,1},{1,0,0},{0,1,0}, {0,0,1},{1,0,0},{0,1,0}};

static void process_tcaf(int nframes,real dt,int nkc,real **tc,rvec *kfac,
			 real rho,real wt,const char *fn_trans,
			 const char *fn_tca,const char *fn_tc,
                         const char *fn_tcf,const char *fn_cub,
			 const char *fn_vk, const output_env_t oenv)
{
  FILE *fp,*fp_vk,*fp_cub=NULL;
  int  nk,ntc;
  real **tcaf,**tcafc=NULL,eta;
  int  i,j,k,kc;
  int  ncorr;
  real fitparms[3],*sig;
  
  nk  = kset_c[nkc];
  ntc = nk*NPK;

  if (fn_trans) {
    fp = xvgropen(fn_trans,"Transverse Current","Time (ps)","TC (nm/ps)",
                  oenv); 
    for(i=0; i<nframes; i++) {
      fprintf(fp,"%g",i*dt);
      for(j=0; j<ntc; j++)
	fprintf(fp," %g",tc[j][i]);
      fprintf(fp,"\n");
    }
    ffclose(fp);
    do_view(oenv,fn_trans,"-nxy");
  }
  
  ncorr = (nframes+1)/2;
  if (ncorr > (int)(5*wt/dt+0.5))
    ncorr = (int)(5*wt/dt+0.5)+1;
  snew(tcaf,nk);
  for(k=0; k<nk; k++)
    snew(tcaf[k],ncorr);
  if (fn_cub) {
     snew(tcafc,nkc);
     for(k=0; k<nkc; k++)
       snew(tcafc[k],ncorr);
  }
  snew(sig,ncorr);
  for(i=0; i<ncorr; i++)
    sig[i]=exp(0.5*i*dt/wt);
  
  low_do_autocorr(fn_tca,oenv,"Transverse Current Autocorrelation Functions",
		  nframes,ntc,ncorr,tc,dt,eacNormal,
		  1,FALSE,FALSE,FALSE,0,0,0,0);
  do_view(oenv,fn_tca,"-nxy");
  
  fp = xvgropen(fn_tc,"Transverse Current Autocorrelation Functions",
		"Time (ps)","TCAF",oenv);
  for(i=0; i<ncorr; i++) {
    kc = 0;
    fprintf(fp,"%g",i*dt);
    for(k=0; k<nk; k++) {
      for(j=0; j<NPK; j++)
	tcaf[k][i] += tc[NPK*k+j][i];
      if (fn_cub) 
	for(j=0; j<NPK; j++)
	  tcafc[kc][i] += tc[NPK*k+j][i];
      if (i == 0)
	fprintf(fp," %g",1.0);
      else {
	tcaf[k][i] /= tcaf[k][0];
	fprintf(fp," %g",tcaf[k][i]);
      }
      if (k+1 == kset_c[kc+1])
	kc++;
    }
    fprintf(fp,"\n");
  }
  ffclose(fp);
  do_view(oenv,fn_tc,"-nxy");
  
  if (fn_cub) {
    fp_cub = xvgropen(fn_cub,"TCAF's and fits", "Time (ps)","TCAF",oenv);
    for(kc=0; kc<nkc; kc++) {
      fprintf(fp_cub,"%g %g\n",0.0,1.0);
      for(i=1; i<ncorr; i++) {
	tcafc[kc][i] /= tcafc[kc][0];
	fprintf(fp_cub,"%g %g\n",i*dt,tcafc[kc][i]);
      }
      fprintf(fp_cub,"&\n");
      tcafc[kc][0] = 1.0;
    }
  }
  
  fp_vk = xvgropen(fn_vk,"Fits","k (nm\\S-1\\N)",
		   "\\8h\\4 (10\\S-3\\N kg m\\S-1\\N s\\S-1\\N)",oenv);
  fprintf(fp_vk,"@    s0 symbol 2\n");
  fprintf(fp_vk,"@    s0 symbol color 1\n");
  fprintf(fp_vk,"@    s0 linestyle 0\n");
  if (fn_cub) {
    fprintf(fp_vk,"@    s1 symbol 3\n");
    fprintf(fp_vk,"@    s1 symbol color 2\n");
  }
  fp = xvgropen(fn_tcf,"TCAF Fits","Time (ps)","",oenv);
  for(k=0; k<nk; k++) {
    tcaf[k][0] = 1.0;
    fitparms[0]  = 1;
    fitparms[1]  = 1;
    do_lmfit(ncorr,tcaf[k],sig,dt,0,0,ncorr*dt,
	     oenv,bDebugMode(),effnVAC,fitparms,0);
    eta = 1000*fitparms[1]*rho/
      (4*fitparms[0]*PICO*norm2(kfac[k])/(NANO*NANO));
    fprintf(stdout,"k %6.3f  tau %6.3f  eta %8.5f 10^-3 kg/(m s)\n",
	    norm(kfac[k]),fitparms[0],eta);
    fprintf(fp_vk,"%6.3f %g\n",norm(kfac[k]),eta);
    for(i=0; i<ncorr; i++)
      fprintf(fp,"%g %g\n",i*dt,fit_function(effnVAC,fitparms,i*dt));
    fprintf(fp,"&\n");
  }
  ffclose(fp);
  do_view(oenv,fn_tcf,"-nxy");

  if (fn_cub) {
    fprintf(stdout,"Averaged over k-vectors:\n");
    fprintf(fp_vk,"&\n");
    for(k=0; k<nkc; k++) {
      tcafc[k][0] = 1.0;
      fitparms[0]  = 1;
      fitparms[1]  = 1;
      do_lmfit(ncorr,tcafc[k],sig,dt,0,0,ncorr*dt,
	       oenv,bDebugMode(),effnVAC,fitparms,0);
      eta = 1000*fitparms[1]*rho/
	(4*fitparms[0]*PICO*norm2(kfac[kset_c[k]])/(NANO*NANO));
      fprintf(stdout,
              "k %6.3f  tau %6.3f  Omega %6.3f  eta %8.5f 10^-3 kg/(m s)\n",
	      norm(kfac[kset_c[k]]),fitparms[0],fitparms[1],eta);
      fprintf(fp_vk,"%6.3f %g\n",norm(kfac[kset_c[k]]),eta);
      for(i=0; i<ncorr; i++)
	fprintf(fp_cub,"%g %g\n",i*dt,fit_function(effnVAC,fitparms,i*dt));
      fprintf(fp_cub,"&\n");
    }
    fprintf(fp_vk,"&\n");
    ffclose(fp_cub);
    do_view(oenv,fn_cub,"-nxy");
  }
  ffclose(fp_vk);
  do_view(oenv,fn_vk,"-nxy");
}


int gmx_tcaf(int argc,char *argv[])
{
  const char *desc[] = {
    "g_tcaf computes tranverse current autocorrelations.",
    "These are used to estimate the shear viscosity eta.",
    "For details see: Palmer, JCP 49 (1994) pp 359-366.[PAR]",
    "Transverse currents are calculated using the",
    "k-vectors (1,0,0) and (2,0,0) each also in the y- and z-direction,",
    "(1,1,0) and (1,-1,0) each also in the 2 other planes (these vectors",
    "are not independent) and (1,1,1) and the 3 other box diagonals (also",
    "not independent). For each k-vector the sine and cosine are used, in",
    "combination with the velocity in 2 perpendicular directions. This gives",
    "a total of 16*2*2=64 transverse currents. One autocorrelation is",
    "calculated fitted for each k-vector, which gives 16 tcaf's. Each of",
    "these tcaf's is fitted to f(t) = exp(-v)(cosh(Wv) + 1/W sinh(Wv)),",
    "v = -t/(2 tau), W = sqrt(1 - 4 tau eta/rho k^2), which gives 16 tau's",
    "and eta's. The fit weights decay with time as exp(-t/wt), the tcaf and",
    "fit are calculated up to time 5*wt.",
    "The eta's should be fitted to 1 - a eta(k) k^2, from which",
    "one can estimate the shear viscosity at k=0.[PAR]",
    "When the box is cubic, one can use the option [TT]-oc[tt], which",
    "averages the tcaf's over all k-vectors with the same length.",
    "This results in more accurate tcaf's.",
    "Both the cubic tcaf's and fits are written to [TT]-oc[tt]",
    "The cubic eta estimates are also written to [TT]-ov[tt].[PAR]",
    "With option [TT]-mol[tt] the transverse current is determined of",
    "molecules instead of atoms. In this case the index group should",
    "consist of molecule numbers instead of atom numbers.[PAR]",
    "The k-dependent viscosities in the [TT]-ov[tt] file should be",
    "fitted to eta(k) = eta0 (1 - a k^2) to obtain the viscosity at",
    "infinite wavelength.[PAR]",
    "NOTE: make sure you write coordinates and velocities often enough.",
    "The initial, non-exponential, part of the autocorrelation function",
    "is very important for obtaining a good fit."
  };
  
  static gmx_bool bMol=FALSE,bK34=FALSE;
  static real wt=5;
  t_pargs pa[] = {
    { "-mol", FALSE, etBOOL, {&bMol},
      "Calculate tcaf of molecules" },
    { "-k34", FALSE, etBOOL, {&bK34},
      "Also use k=(3,0,0) and k=(4,0,0)" },
    { "-wt", FALSE, etREAL, {&wt},
      "Exponential decay time for the TCAF fit weights" }
  };

  t_topology top;
  int        ePBC;
  t_trxframe fr;
  matrix     box;
  gmx_bool       bTPS,bTop; /* ,bCubic; */
  int        gnx;
  atom_id    *index,*atndx=NULL,at;
  char       *grpname;
  char       title[256];
  real       t0,t1,dt,m,mtot,sysmass,rho,sx,cx;
  t_trxstatus *status;
  int        nframes,n_alloc,i,j,k,d;
  rvec       mv_mol,cm_mol,kfac[NK];
  int        nkc,nk,ntc;
  real       **c1,**tc;
  output_env_t oenv;
  
#define NHISTO 360
    
  t_filenm  fnm[] = {
    { efTRN, "-f",    NULL,      ffREAD  },
    { efTPS, NULL,    NULL,      ffOPTRD }, 
    { efNDX, NULL,    NULL,      ffOPTRD },
    { efXVG, "-ot",  "transcur", ffOPTWR },
    { efXVG, "-oa",  "tcaf_all", ffWRITE },
    { efXVG, "-o",   "tcaf",     ffWRITE },
    { efXVG, "-of",  "tcaf_fit", ffWRITE },
    { efXVG, "-oc",  "tcaf_cub", ffOPTWR },
    { efXVG, "-ov",  "visc_k",   ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
    
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL,&oenv);

  bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,NULL,NULL,box,
		     TRUE);
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

  if (bMol) {
    if (!bTop)
      gmx_fatal(FARGS,"Need a topology to determine the molecules");
    atndx = top.mols.index;
  }

  if (bK34)
    nkc = NKC;
  else
    nkc = NKC0;
  nk  = kset_c[nkc];
  ntc = nk*NPK; 

  sprintf(title,"Velocity Autocorrelation Function for %s",grpname);
  
  sysmass = 0;
  for(i=0; i<nk; i++) {
    if (iprod(v0[i],v1[i]) != 0)
      gmx_fatal(FARGS,"DEATH HORROR: vectors not orthogonal");
    if (iprod(v0[i],v2[i]) != 0)
      gmx_fatal(FARGS,"DEATH HORROR: vectors not orthogonal");
    if (iprod(v1[i],v2[i]) != 0)
	gmx_fatal(FARGS,"DEATH HORROR: vectors not orthogonal");
    unitv(v1[i],v1[i]);
    unitv(v2[i],v2[i]);
  }
  snew(tc,ntc);
    for(i=0; i<top.atoms.nr; i++)
      sysmass += top.atoms.atom[i].m;

  read_first_frame(oenv,&status,ftp2fn(efTRN,NFILE,fnm),&fr,
		   TRX_NEED_X | TRX_NEED_V);
  t0=fr.time;
  
  n_alloc=0;
  nframes=0;
  rho=0;
  /* bCubic = TRUE; */
  do {
    /*
    bCubic = bCubic && !TRICLINIC(fr.box) &&
    fabs(fr.box[XX][XX]-fr.box[YY][YY]) < 0.001*fr.box[XX][XX] &&
      fabs(fr.box[XX][XX]-fr.box[ZZ][ZZ]) < 0.001*fr.box[XX][XX];
      */

    if (nframes >= n_alloc) {
      n_alloc+=100;
      for (i=0; i<ntc; i++)
	srenew(tc[i],n_alloc);
    }

    rho += 1/det(fr.box);
    for(k=0; k<nk; k++)
      for(d=0; d<DIM; d++)
	kfac[k][d] = 2*M_PI*v0[k][d]/fr.box[d][d];
    for (i=0; i<ntc; i++)
      tc[i][nframes] = 0;

    for(i=0; i<gnx; i++) {
      if (bMol) {
	clear_rvec(mv_mol);
	clear_rvec(cm_mol);
	mtot = 0;
	for(j=0; j<atndx[index[i]+1] - atndx[index[i]]; j++) {
	  at = atndx[index[i]] + j;
	  m  = top.atoms.atom[at].m;
	  mv_mol[XX] += m*fr.v[at][XX];
	  mv_mol[YY] += m*fr.v[at][YY];
	  mv_mol[ZZ] += m*fr.v[at][ZZ];
	  cm_mol[XX] += m*fr.x[at][XX];
	  cm_mol[YY] += m*fr.x[at][YY];
	  cm_mol[ZZ] += m*fr.x[at][ZZ];
	  mtot += m;
	}
	svmul(1.0/mtot,cm_mol,cm_mol);
      } else
	svmul(top.atoms.atom[index[i]].m,fr.v[index[i]],mv_mol);
      
      if (!bMol)
	copy_rvec(fr.x[index[i]],cm_mol);
      j=0;
      for(k=0; k<nk; k++) {
	sx = sin(iprod(kfac[k],cm_mol));
	cx = cos(iprod(kfac[k],cm_mol));
	tc[j][nframes] += sx*iprod(v1[k],mv_mol);
	j++;
	tc[j][nframes] += cx*iprod(v1[k],mv_mol);
	j++;
	tc[j][nframes] += sx*iprod(v2[k],mv_mol);
	j++;
	tc[j][nframes] += cx*iprod(v2[k],mv_mol);
	j++;
      }
    }
    
    t1=fr.time;
    nframes ++;
  } while (read_next_frame(oenv,status,&fr));
  close_trj(status);

  dt = (t1-t0)/(nframes-1);

  rho *= sysmass/nframes*AMU/(NANO*NANO*NANO);
  fprintf(stdout,"Density = %g (kg/m^3)\n",rho);
  process_tcaf(nframes,dt,nkc,tc,kfac,rho,wt,
	       opt2fn_null("-ot",NFILE,fnm),
	       opt2fn("-oa",NFILE,fnm),opt2fn("-o",NFILE,fnm),
	       opt2fn("-of",NFILE,fnm),opt2fn_null("-oc",NFILE,fnm),
	       opt2fn("-ov",NFILE,fnm),oenv);
  
  thanx(stderr);
  
  return 0;
}
