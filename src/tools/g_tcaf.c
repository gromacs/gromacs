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

#define NK  24
#define NPK 4

#define NKC  6
#define NKC0 4
int kset_c[NKC+1] = { 0, 3, 9, 13, 16, 19, NK };

rvec v0[NK]={{1,0,0},{0,1,0},{0,0,1}, {1, 1,0},{1,-1,0},{1,0,1}, {1,0,-1},{0,1, 1},{0,1,-1}, {1, 1, 1},{1, 1,-1},{1,-1, 1},{-1,1, 1}, {2,0,0},{0,2,0},{0,0,2}, {3,0,0},{0,3,0},{0,0,3}, {4,0,0},{0,4,0},{0,0,4}};
rvec v1[NK]={{0,1,0},{0,0,1},{1,0,0}, {0, 0,1},{0, 0,1},{0,1,0}, {0,1, 0},{1,0, 0},{1,0, 0}, {1,-1, 0},{1,-1, 0},{1, 0,-1},{ 0,1,-1}, {0,1,0},{0,0,1},{1,0,0}, {0,1,0},{0,0,1},{1,0,0}, {0,1,0},{0,0,1},{1,0,0}};
rvec v2[NK]={{0,0,1},{1,0,0},{0,1,0}, {1,-1,0},{1, 1,0},{1,0,-1},{1,0, 1},{0,1,-1},{0,1, 1}, {1, 1,-2},{1, 1, 2},{1, 2, 1},{ 2,1, 1}, {0,0,1},{1,0,0},{0,1,0}, {0,0,1},{1,0,0},{0,1,0}, {0,0,1},{1,0,0},{0,1,0}};

static void process_tcaf(int nframes,real dt,int nkc,real **tc,rvec *kfac,
			 real rho,real endfit,char *fn_trans,
			 char *fn_tca,char *fn_tc,char *fn_tcf,char *fn_cub,
			 char *fn_vk)
{
  FILE *fp,*fp_vk,*fp_cub;
  int  nk,ntc;
  real **tcaf,**tcafc,eta;
  int  i,j,k,kc;
  int  ncorr;
  real fitparms[2],*sig,*fit;
  
  nk  = kset_c[nkc];
  ntc = nk*NPK;

  if (fn_trans) {
    fp = xvgropen(fn_trans,"Transverse Current","Time (ps)","TC (nm/ps)"); 
    for(i=0; i<nframes; i++) {
      fprintf(fp,"%g",i*dt);
      for(j=0; j<ntc; j++)
	fprintf(fp," %g",tc[j][i]);
      fprintf(fp,"\n");
    }
    fclose(fp);
    do_view(fn_trans,NULL);
  }
  
  ncorr = (int)(endfit/dt+0.5)+1;
  if (ncorr > (nframes+1)/2)
    ncorr = (nframes+1)/2;
  snew(tcaf,nk);
  for(k=0; k<nk; k++)
    snew(tcaf[k],ncorr);
  if (fn_cub) {
     snew(tcafc,nkc);
     for(k=0; k<nkc; k++)
       snew(tcafc[k],ncorr);
  }
  snew(sig,ncorr);
  snew(fit,ncorr);
  sig[0] = 1;
  for(i=1; i<ncorr; i++)
    sig[i]=sqrt(i);
  
  low_do_autocorr(fn_tca,"Transverse Current Autocorrelation Functions",
		  nframes,ntc,-1,tc,dt,eacNormal,
		  1,FALSE,TRUE,FALSE,FALSE,0,0,0,0);
  do_view(fn_tca,NULL);
  
  fp = xvgropen(fn_tc,"Transverse Current Autocorrelation Functions",
		"Time (ps)","TCAF");
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
  fclose(fp);
  do_view(fn_tc,NULL);
  
  if (fn_cub) {
    fp_cub = xvgropen(fn_cub,"TCAF's and fits",
		      "Time (ps)","TCAF");
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
		   "eta (10\\S-3\\N kg m\\S-1\\N s\\S-1\\N)");
  fprintf(fp_vk,"@    s0 symbol 2\n");
  fprintf(fp_vk,"@    s0 symbol color 1\n");
  fprintf(fp_vk,"@    s0 linestyle 0\n");
  if (fn_cub) {
    fprintf(fp_vk,"@    s1 symbol 3\n");
    fprintf(fp_vk,"@    s1 symbol color 2\n");
  }
  fp = xvgropen(fn_tcf,"TCAF Fits","Time (ps)","");
  for(k=0; k<nk; k++) {
    tcaf[k][0] = 1.0;
    fitparms[0]  = 1;
    fitparms[1]  = 1;
    do_lmfit(ncorr,tcaf[k],sig,dt,0,0,endfit,
	     bDebugMode(),effnVAC,fitparms,0);
    eta = 1000*fitparms[1]*rho/
      (4*fitparms[0]*PICO*norm2(kfac[k])/(NANO*NANO));
    fprintf(stdout,"k %6.3f  tau %6.3f  eta %8.5f 10^-3 kg/(m s)\n",
	    norm(kfac[k]),fitparms[0],eta);
    fprintf(fp_vk,"%6.3f %g\n",norm(kfac[k]),eta);
    for(i=0; i<ncorr; i++)
      fprintf(fp,"%g %g\n",i*dt,fit_function(effnVAC,fitparms,i*dt));
    fprintf(fp,"&\n");
  }
  fclose(fp);
  do_view(fn_tcf,NULL);

  if (fn_cub) {
    fprintf(stdout,"Averaged over k-vectors:\n");
    fprintf(fp_vk,"&\n");
    for(k=0; k<nkc; k++) {
      tcafc[k][0] = 1.0;
      fitparms[0]  = 1;
      fitparms[1]  = 1;
      do_lmfit(ncorr,tcafc[k],sig,dt,0,0,endfit,
	       bDebugMode(),effnVAC,fitparms,0);
      eta = 1000*fitparms[1]*rho/
	(4*fitparms[0]*PICO*norm2(kfac[kset_c[k]])/(NANO*NANO));
      fprintf(stdout,"k %6.3f  tau %6.3f  eta %8.5f 10^-3 kg/(m s)\n",
	      norm(kfac[kset_c[k]]),fitparms[0],eta);
      fprintf(fp_vk,"%6.3f %g\n",norm(kfac[kset_c[k]]),eta);
      for(i=0; i<ncorr; i++)
	fprintf(fp_cub,"%g %g\n",i*dt,fit_function(effnVAC,fitparms,i*dt));
      fprintf(fp_cub,"&\n");
    }
    fprintf(fp_vk,"&\n");
    fclose(fp_cub);
    do_view(fn_cub,NULL);
  }
  fclose(fp_vk);
  do_view(fn_vk,NULL);

  
}


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_tcaf computes tranverse current autocorrelations.",
    "These are used to estimate the shear viscosity eta.",
    "Transverse currents are calculated using the",
    "k-vectors (1,0,0) and (2,0,0) each also in the y- and z-direction,",
    "(1,1,0) and (1,-1,0) each also in the 2 other plains (these vectors",
    "are not independent) and (1,1,1) and the 3 other box diagonals (also",
    "not independent). For each k-vector the sine and cosine are used, in",
    "combination with the velocity in 2 perpendicular directions. This gives",
    "a total of 16*2*2=64 transverse currents. One autocorrelation is",
    "calculated fitted for each k-vector, which gives 16 tcaf's. Each of",
    "these tcaf's is fitted to f(t) = exp(-v)(cosh(wv) + 1/w sinh(wv)),",
    "v = -t/(2 tau), w = sqrt(1 - 4 tau eta/rho k^2), which gives 16 tau's",
    "and eta's. The fit weights decay with time as 1/t.",
    "The eta's should be fitted to 1 - a eta(k) k^2, from which",
    "one can estimate the shear viscosity at k=0.[PAR]",
    "When the box is cubic, one can use the option [TT]-oc[tt], which",
    "averages the tcaf's over all k-vectors with the same length.",
    "This results in more accurate tcaf's.",
    "Both the cubic tcaf's and fits are written to [TT]-oc[tt]",
    "The cubic eta estimates are also written to [TT]-vk[tt].[PAR]",
    "With option [TT]-mol[tt] the transverse current is determined of",
    "molecules instead of atoms. In this case the index group should",
    "consist of molecule numbers instead of atom numbers.",
  };
  
  static bool bMol=FALSE,bK34=FALSE;
  static real endfit=50;
  t_pargs pa[] = {
    { "-mol", FALSE, etBOOL, {&bMol},
      "Calculate tcaf of molecules" },
    { "-k34", FALSE, etBOOL, {&bK34},
      "Also use k=(3,0,0) and k=(4,0,0)" },
    { "-endfit", FALSE, etREAL, {&endfit},
      "Time where to end the tcaf plots and fits" }
  };

  t_topology top;
  t_trxframe fr;
  matrix     box;
  bool       bTPS,bTop; /* ,bCubic; */
  int        gnx;
  atom_id    *index,*a,*atndx,at;
  char       *grpname;
  char       title[256];
  real       t0,t1,dt,m,mtot,sysmass,rho,sx,cx;
  int        status,nframes,n_alloc,i,j,k,d,tel3;
  rvec       mv_mol,cm_mol,kfac[NK];
  int        nkc,nk,ntc;
  real       **c1,**tc;
  
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

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL);

  bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,NULL,NULL,box,TRUE);
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

  if (bMol) {
    if (!bTop)
      fatal_error(0,"Need a topology to determine the molecules");
    a     = top.blocks[ebMOLS].a;
    atndx = top.blocks[ebMOLS].index;
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
      fatal_error(0,"DEATH HORROR: vectors not orthogonal");
    if (iprod(v0[i],v2[i]) != 0)
      fatal_error(0,"DEATH HORROR: vectors not orthogonal");
    if (iprod(v1[i],v2[i]) != 0)
	fatal_error(0,"DEATH HORROR: vectors not orthogonal");
    unitv(v1[i],v1[i]);
    unitv(v2[i],v2[i]);
  }
  snew(tc,ntc);
    for(i=0; i<top.atoms.nr; i++)
      sysmass += top.atoms.atom[i].m;

  read_first_frame(&status,ftp2fn(efTRN,NFILE,fnm),&fr,
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
    tel3=3*nframes;

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
	  at = a[atndx[index[i]]+j];
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
  } while (read_next_frame(status,&fr));
  close_trj(status);

  dt = (t1-t0)/(nframes-1);

  rho *= sysmass/nframes*AMU/(NANO*NANO*NANO);
  fprintf(stdout,"Density = %g (kg/m^3)\n",rho);
  process_tcaf(nframes,dt,nkc,tc,kfac,rho,endfit,
	       opt2fn_null("-ot",NFILE,fnm),
	       opt2fn("-oa",NFILE,fnm),opt2fn("-o",NFILE,fnm),
	       opt2fn("-of",NFILE,fnm),opt2fn_null("-oc",NFILE,fnm),
	       opt2fn("-ov",NFILE,fnm));
  
  thanx(stderr);
  
  return 0;
}
