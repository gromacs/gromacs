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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_g_velacc_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "confio.h"
#include "copyrite.h"
#include "fatal.h"
#include "futil.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "rdgroup.h"
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

#define NK  16
#define NPK 4
#define NTC (NK*NPK)

int kset_r[NK+1] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 };

#define NKC 4
int kset_c[NKC+1] = { 0, 3, 9, 13, NK };

rvec v0[NK]={{1,0,0},{0,1,0},{0,0,1}, {1, 1,0},{1,-1,0},{1,0,1}, {1,0,-1},{0,1, 1},{0,1,-1}, {1, 1, 1},{1, 1,-1},{1,-1, 1},{-1,1, 1}, {2,0,0},{0,2,0},{0,0,2}};
rvec v1[NK]={{0,1,0},{0,0,1},{1,0,0}, {0, 0,1},{0, 0,1},{0,1,0}, {0,1, 0},{1,0, 0},{1,0, 0}, {1,-1, 0},{1,-1, 0},{1, 0,-1},{ 0,1,-1}, {0,1,0},{0,0,1},{1,0,0}};
rvec v2[NK]={{0,0,1},{1,0,0},{0,1,0}, {1,-1,0},{1, 1,0},{1,0,-1},{1,0, 1},{0,1,-1},{0,1, 1}, {1, 1,-2},{1, 1, 2},{1, 2, 1},{ 2,1, 1}, {0,0,1},{1,0,0},{0,1,0}};

static void process_tcaf(int teller,real dt,real **tc,rvec *kfac,
			 real rho,real wt,bool bCubic,
			 char *fn_tca,char *fn_tc,char *fn_tcf,char *fn_vk)
{
  FILE *fp;
  FILE *fp2;
  real **tcaf,eta,*eta_av;
  int  i,j,k,k0,nk,*kset;
  int  ncorr;
  real fitparms[2],*sig,*fit;

  if (bCubic) {
    nk = NKC;
    kset = kset_c;
  } else {
    nk = NK;
    kset = kset_r;
  }

  if (bDebugMode) {
    fp = xvgropen("transcur.xvg", 
		  "Transverse Current","Time (ps)","TC (nm/ps)"); 
    for(i=0; i<teller; i++) {
      fprintf(fp,"%g",i*dt);
      for(j=0; j<NTC; j++)
	fprintf(fp," %g",tc[j][i]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  
  ncorr = (teller+1)/2;
  if (ncorr > (int)(4*wt/dt+0.5))
    ncorr = (int)(4*wt/dt+0.5)+1;
  snew(tcaf,nk);
  for(k0=0; k0<nk; k0++)
    snew(tcaf[k0],ncorr);
  snew(sig,ncorr);
  snew(fit,ncorr);
  for(i=0; i<ncorr; i++)
    sig[i]=exp(i*dt/wt);
  
  low_do_autocorr(fn_tca,"Transverse Current Autocorrelation Functions",
		  teller,NTC,-1,tc,dt,eacNormal,
		  1,FALSE,TRUE,FALSE,FALSE,0,0,0,0);
  do_view(fn_tca,NULL);
  
  fp = xvgropen(fn_tc,"Transverse Current Autocorrelation Functions",
		"Time (ps)","TCAF");
  for(i=0; i<ncorr; i++) {
    k0 = 0;
    fprintf(fp,"%g",i*dt);
    for(k=0; k<NK; k++) {
      for(j=0; j<NPK; j++)
	tcaf[k0][i] += tc[NPK*k+j][i];
      if (k+1 == kset[k0+1]) {
	if (i == 0)
	  fprintf(fp," %g",1.0);
	else {
	  tcaf[k0][i] /= tcaf[k0][0];
	  fprintf(fp," %g",tcaf[k0][i]);
	}
	k0++;
      }
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  do_view(fn_tc,NULL);
  
  fp = xvgropen(fn_vk,"Fits","k (nm\\S-1\\N)",
		"eta (kg m\\S-1\\N s\\S-1\\N)");
  fp2 = xvgropen(fn_tcf,"TCAF Fits","Time (ps)","");
  k0 = 0;
  for(k=0; k<nk; k++) {
    tcaf[k][0] = 1.0;
    fitparms[0]  = 1;
    fitparms[1]  = 1;
    do_lmfit(ncorr,tcaf[k],sig,dt,0,0,ncorr*dt,
	     bDebugMode(),effnVAC,fitparms,fit,NULL);
    eta = 1000*fitparms[1]*rho/
      (4*fitparms[0]*PICO*norm2(kfac[kset[k]])/(NANO*NANO));
    fprintf(stdout,"k %6.3f  tau %6.3f  eta %8.5f 10^-3 kg/(m s)\n",
	    norm(kfac[kset[k]]),fitparms[0],eta);
    fprintf(fp,"%6.3f %g\n",norm(kfac[kset[k]]),eta);
    for(i=0; i<ncorr; i++)
      fprintf(fp2,"%g %g\n",i*dt,fit[i]);
    fprintf(fp2,"&\n");
  }
  fclose(fp2);
  fclose(fp);
  do_view(fn_tcf,NULL);
  do_view(fn_vk,NULL);
}


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_velacc computes the velocity autocorrelation function.",
    "When the [TT]-s[tt] option is used, the momentum autocorrelation",
    "function is calculated.[PAR]",
    "With option [TT]-mol[tt] the momentum autocorrelation function of",
    "molecules is calculated. In this case the index group should consist",
    "of molecule numbers instead of atom numbers.[PAR]",
    "Setting any of the last 4 output file options will turn on the",
    "tranverse current autocorrelation calculation. This always produces",
    "all the 4 output files. Transverse currents are calculated using the",
    "k-vectors (1,0,0) and (2,0,0) each also in the y- and z-direction,",
    "(1,1,0) and (1,-1,0) each also in the 2 other plains (these vectors",
    "are not independent) and (1,1,1) and the 3 other box diagonals (also",
    "not independent). For each k-vector the sine and cosine are used, in",
    "combination with the velocity in 2 perpendicular directions. This gives",
    "a total of 16*2*2=64 transverse currents. One autocorrelation is",
    "calculated fitted for each k-vector, which gives 16 tcaf's. Each of",
    "these tcaf's is fitted to f(t) = exp(-v)(cosh(wv) + 1/w sinh(wv)),",
    "v = -t/(2 tau), w = sqrt(1 - 4 tau eta/rho k^2), which gives 16 tau's",
    "and eta's. The eta's should be fitted to 1 - a eta(k) k^2, from which",
    "one can estimate the shear viscosity at k=0."
  };
  
  static bool bMol=FALSE,bCubic=FALSE;
  static real wt=5;
  t_pargs pa[] = {
    { "-mol", FALSE, etBOOL, {&bMol},
      "Calculate vac of molecules" },
    { "-cubic", FALSE, etBOOL, {&bCubic},
      "Assume a cubic box for tcaf calculation" },
    { "-wt", FALSE, etREAL, {&wt},
      "Exponential decay time for the TCAF fit weights" }
  };

  t_topology top;
  t_trxframe fr;
  matrix     box;
  bool       bTPS,bTop,bVAC,bTCAF; /* ,bCubic; */
  int        gnx;
  atom_id    *index,*a,*atndx,at;
  char       *grpname;
  char       title[256];
  real       t0,t1,dt,m,mtot,sysmass,rho,sx,cx;
  int        status,teller,n_alloc,i,j,k,d,tel3;
  rvec       mv_mol,cm_mol,kfac[NK];
  real       **c1,**tc;
  
#define NHISTO 360
    
  t_filenm  fnm[] = {
    { efTRN, "-f",    NULL,   ffREAD  },
    { efTPS, NULL,    NULL,   ffOPTRD }, 
    { efNDX, NULL,    NULL,   ffOPTRD },
    { efXVG, "-o",    "vac",  ffOPTWR },
    { efXVG, "-tc",   "tcaf", ffOPTWR },
    { efXVG, "-tca",  "tcaf_all", ffOPTWR },
    { efXVG, "-tcf",  "tcaf_fit", ffOPTWR },
    { efXVG, "-vk",   "visc_k", ffOPTWR }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);

  bVAC  = opt2bSet("-o",NFILE,fnm);
  bTCAF = opt2bSet("-tc",NFILE,fnm) || opt2bSet("-tca",NFILE,fnm)
    || opt2bSet("-tcf",NFILE,fnm) || opt2bSet("-vk",NFILE,fnm);

  if (!bVAC && !bTCAF) {
    fprintf(stderr,"Nothing to do\n");
    exit(0);
  }

  bTPS = bMol || ftp2bSet(efTPS,NFILE,fnm) || !ftp2bSet(efNDX,NFILE,fnm)
    || bTCAF;

  if (bTPS) {
    bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,NULL,NULL,box,TRUE);
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  } else
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

  if (bMol) {
    if (!bTop)
      fatal_error(0,"Need a topology to determine the molecules");
    a     = top.blocks[ebMOLS].a;
    atndx = top.blocks[ebMOLS].index;
  }
  
  sprintf(title,"Velocity Autocorrelation Function for %s",grpname);
  
  if (bVAC)
    snew(c1,gnx);

  sysmass = 0;
  if (bTCAF) {
    for(i=0; i<NK; i++) {
      if (iprod(v0[i],v1[i]) != 0)
	fatal_error(0,"DEATH HORROR: vectors not orthogonal");
      if (iprod(v0[i],v2[i]) != 0)
	fatal_error(0,"DEATH HORROR: vectors not orthogonal");
      if (iprod(v1[i],v2[i]) != 0)
	fatal_error(0,"DEATH HORROR: vectors not orthogonal");
      unitv(v1[i],v1[i]);
      unitv(v2[i],v2[i]);
    }
    snew(tc,NTC);
    for(i=0; i<top.atoms.nr; i++)
      sysmass += top.atoms.atom[i].m;
  }

  read_first_frame(&status,ftp2fn(efTRN,NFILE,fnm),&fr,
		   bTCAF ? TRX_READ_X | TRX_NEED_V  : TRX_NEED_V);
  t0=fr.time;
      
  n_alloc=0;
  teller=0;
  rho=0;
  /* bCubic = TRUE; */
  do {
    /*
    bCubic = bCubic && !TRICLINIC(fr.box) &&
      fabs(fr.box[XX][XX]-fr.box[YY][YY]) < 0.001*fr.box[XX][XX] &&
      fabs(fr.box[XX][XX]-fr.box[ZZ][ZZ]) < 0.001*fr.box[XX][XX];
      */

    if (teller >= n_alloc) {
      n_alloc+=100;
      if (bVAC)
	for(i=0; i<gnx; i++)
	  srenew(c1[i],DIM*n_alloc);
      if (bTCAF)
	for (i=0; i<NTC; i++)
	  srenew(tc[i],n_alloc);
    }
    tel3=3*teller;

    if (bTCAF) {
      rho += 1/det(fr.box);
      for(k=0; k<NK; k++)
	for(d=0; d<DIM; d++)
	  kfac[k][d] = 2*M_PI*v0[k][d]/fr.box[d][d];
      for (i=0; i<NTC; i++)
	tc[i][teller] = 0;
    }

    for(i=0; i<gnx; i++) {
      if (bMol) {
	clear_rvec(mv_mol);
	clear_rvec(cm_mol);
	mtot = 0;
	for(j=0; j<atndx[index[i]+1] - atndx[index[i]]; j++) {
	  at = a[atndx[index[i]]+j];
	  m  = top.atoms.atom[at].m;
	  mv_mol[XX] += m*fr.v[at][XX];
	  mv_mol[YY] += m*fr.v[at][YY];
	  mv_mol[ZZ] += m*fr.v[at][ZZ];
	  if (fr.bX) {
	    mtot += m;
	    cm_mol[XX] += m*fr.x[at][XX];
	    cm_mol[YY] += m*fr.x[at][YY];
	    cm_mol[ZZ] += m*fr.x[at][ZZ];
	  }
	}
	if (fr.bX)
	  svmul(1.0/mtot,cm_mol,cm_mol);
      } else {
	if (bTPS)
	  m = top.atoms.atom[index[i]].m;
	else
	  m = 1;
	svmul(m,fr.v[index[i]],mv_mol);
      }
      if (bVAC) {
	c1[i][tel3+XX]=mv_mol[XX];
	c1[i][tel3+YY]=mv_mol[YY];
	c1[i][tel3+ZZ]=mv_mol[ZZ];
      }
      if (bTCAF && fr.bX) {
	if (!bMol)
	  copy_rvec(fr.x[index[i]],cm_mol);
	j=0;
	for(k=0; k<NK; k++) {
	  sx = sin(iprod(kfac[k],cm_mol));
	  cx = cos(iprod(kfac[k],cm_mol));
	  tc[j][teller] += sx*iprod(v1[k],mv_mol);
	  j++;
	  tc[j][teller] += cx*iprod(v1[k],mv_mol);
	  j++;
	  tc[j][teller] += sx*iprod(v2[k],mv_mol);
	  j++;
	  tc[j][teller] += cx*iprod(v2[k],mv_mol);
	  j++;
	}
      }
    }

    t1=fr.time;
    teller ++;
  } while (read_next_frame(status,&fr));
  close_trj(status);

  dt = (t1-t0)/(teller-1);

  if (bVAC) {
    do_autocorr(opt2fn("-o",NFILE,fnm),"Velocity Autocorrelation Function",
		teller,gnx,c1,dt,eacVector,TRUE);
    do_view(opt2fn("-o",NFILE,fnm),"-nxy");
  }

  if (bTCAF) {
    rho *= sysmass/teller*AMU/(NANO*NANO*NANO);
    fprintf(stdout,"Density = %g (kg/m^3)\n",rho);
    process_tcaf(teller,dt,tc,kfac,rho,wt,bCubic,
		 opt2fn("-tca",NFILE,fnm),opt2fn("-tc",NFILE,fnm),
		 opt2fn("-tcf",NFILE,fnm),opt2fn("-vk",NFILE,fnm));
  }
  
  thanx(stderr);
  
  return 0;
}
