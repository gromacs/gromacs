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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_nrnb_c = "$Id$";

#include <sysstuff.h>
#include "fatal.h"
#include "vveclib.h"
#include "names.h"
#include "macros.h"
#include "nrnb.h"
#include "main.h"
#include "smalloc.h"
#include "copyrite.h"

typedef struct {
  char *name;
  int  flop;
  int  cycles;
} t_nrnb_data;

static t_nrnb_data nbdata[eNRNB] = {
  { "LJ+Coulomb",    38, 118 },
  { "Coulomb",       27,  83 },
  { "LJC-RF",        43, 130 },
  { "Coul-RF",       32,  90 },
  { "FreeEnergyLJC", 64, 110 },
  { "Table-Coul",    39,  25 },
  { "Table-LJC",     47,  30 },
  { "Calc Weights",  36,  80 },
  { "Spread Q",       6,  21 },
  { "Gather F",      23,  46 },
  { "3D-FFT",         8,  20 },
  { "Convolution",    4,  10 },
  { "Long-Range",    27,  83 },
  { "Buckingham",    60, 150 },
  { "Buck.+RF",      65, 180 },
  { "NS-Pairs",      21, 144 },
  { "Reset In Box",   9,  53 },
  { "Shift-X",        6,  63 },
  { "CG-CoM",        29,  50 },
  { "Sum Forces",     1,   3 },
  { "Bonds",         43,  90 },
  { "Angles",       163, 325 },
  { "Propers",      229, 460 },
  { "Impropers",    208, 420 },
  { "RB-Dihedrals", 247, 500 },
  { "Dist. Restr.", 200, 500 },
  { "Pos. Restr.",   50, 100 },
  { "Morse Potent.",  0,   0 },
  { "Virial",        18,  43 },
  { "Update",        31, 273 },
  { "Stop-CM",       10,  62 },
  { "P-Coupling",    24,  58 },
  { "Calc-Ekin",     27,  18 },
  { "Shake",         30,  99 },
  { "Shake-V",       15,  25 },
  { "Shake-Init",    10,  52 },
  { "Shake-Vir",     18,  43 },
  { "Settle",       323, 700 },
  { "PShake-InitLD", 59, 177 },    
  { "PShake-InitMD", 65, 195 },   
  { "PShake",         7,  25 },
  { "Dummy1",        17,  70 },
  { "Dummy2",        28,  80 },
  { "Dummy2fd",      73, 250 }, /* 250 is a wild guess */
  { "Dummy2fad",    131, 250 }, /* 250 is a wild guess */
  { "Dummy3",        87, 300 },
};

void init_nrnb(t_nrnb *nrnb)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    nrnb->n[i]=0.0;
}

void cp_nrnb(t_nrnb *dest, t_nrnb *src)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    dest->n[i]=src->n[i];
}

void add_nrnb(t_nrnb *dest, t_nrnb *s1, t_nrnb *s2)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    dest->n[i]=s1->n[i]+s2->n[i];
}

void print_nrnb(FILE *out, t_nrnb *nrnb)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    fprintf(out,"%14s  %10.0f.\n",nbdata[i].name,nrnb->n[i]);
}

void _inc_nrnb(t_nrnb *nrnb,int enr,int inc,char *file,int line)
{
  nrnb->n[enr]+=inc;
#ifdef DEBUG
  fprintf(stdlog,"nrnb %15s(%2d) incremented with %8d from file %s line %d\n",
	  nbdata[enr].name,enr,inc,file,line);
#endif
}

void print_perf(FILE *out,double cputime,double realtime,t_nrnb *nrnb,int nprocs)
{
  int    nbfs_ind[] = { eNR_LJC, eNR_QQ, eNR_LJCRF, eNR_QQRF, 
			eNR_FREE, eNR_COULTAB, eNR_TAB, eNR_LR, 
			eNR_BHAM, eNR_BHAMRF };
  int    i;
  double nbfs,mni,frac,tfrac;
  double mflop,dni,tni,tcyc;
  
  if (cputime == 0.0) {
    fprintf(out,"cputime = 0! Infinite Giga flopses! \n");
    return;
  }
  nbfs=0.0;
  for(i=0; (i<asize(nbfs_ind)); i++) {
    nbfs+=1e-6*nrnb->n[nbfs_ind[i]];
  }
  tcyc=0;
  for(i=0; (i<eNRNB); i++) {
    tcyc+=1e-6*nrnb->n[i]*nbdata[i].cycles;
  }
  if (tcyc == 0) {
    fprintf(out,"No MEGA Flopsen this time\n");
    return;
  }
  fprintf(out,"\tM E G A - F L O P S   A C C O U N T I N G\n\n");
  fprintf(out,"%15s  %12s  %12s  %12s  %8s\n",
	  "Computing:","M-Number","M-Flops","M-Cycles","% Time");
  mflop=0.0;
  tni=0.0;
  tfrac=0.0;
  for(i=0; (i<eNRNB); i++) {
    mni    = 1e-6*nrnb->n[i];
    mflop += mni*nbdata[i].flop;
    dni    = mni*nbdata[i].cycles;
    frac   = 100.0*mni*nbdata[i].cycles/tcyc;
    tfrac += frac;
    tni   += dni;
    if (mni != 0)
      fprintf(out,"%15s  %12.5f  %12.5f  %12.5f  %6.1f\n",
	      nbdata[i].name,mni,mni*nbdata[i].flop,dni,frac);
  }
  fprintf(out,"%15s  %12s  %12.5f  %12.5f  %6.1f\n\n",
	  "Total","",mflop,tni,tfrac);
  fprintf(out,"CPU time:    %10.3f s.\n",cputime);
  fprintf(out,"Real time:   %10.3f s. [%.1f%%]\n",
	  realtime,100.0*cputime/realtime);
  if (cputime > 60) {
    fprintf(out,"             ");
    pr_difftime(out,cputime);
  }
  fprintf(out,"Performance: %10.3f Mnbf/s\n",nbfs/cputime);
  fprintf(out,"             %10.3f MFlops\n",mflop/cputime);
}

int cost_nrnb(int enr)
{
  return nbdata[enr].cycles;
}

char *nrnb_str(int enr)
{
  return nbdata[enr].name;
}

static int    force_index[]={ 
  eNR_LJC,    eNR_QQ,     eNR_BHAM, 
  eNR_BONDS,  eNR_ANGLES,  eNR_PROPER, eNR_IMPROPER, 
  eNR_RB,     eNR_DISRES,  eNR_POSRES,
  eNR_NS
};
#define NFORCE_INDEX asize(force_index)

static int    shake_index[]={ 
  eNR_SHAKE,  eNR_SHAKE_RIJ, eNR_SETTLE, eNR_UPDATE, eNR_PCOUPL,
  eNR_SHAKE_VIR, eNR_SHAKE_V
};
#define NSHAKE_INDEX asize(shake_index)

static double pr_av(FILE *log,int nprocs,double fav,double ftot[],char *title)
{
  int    i,perc;
  double dperc,unb;
  
  unb=0;
  if (fav > 0) {
    fav/=nprocs;
    fprintf(log,"\n%15s:",title);
    for(i=0; (i<nprocs); i++) {
      dperc=(100.0*ftot[i])/fav;
      unb=max(unb,dperc);
      perc=dperc;
      fprintf(log,"%3d ",perc);
    }
    if (unb > 0) {
      perc=10000.0/unb;
      fprintf(log,"%6d%%\n\n",perc);
    }
    else
      fprintf(log,"\n\n");
  }
  return unb;
}

void pr_load(FILE *log,int nprocs,t_nrnb nrnb[])
{
  int    i,j,perc;
  double dperc,unb,uf,us;
  double *ftot,fav;
  double *stot,sav;
  t_nrnb *av;

  snew(av,1);
  snew(ftot,nprocs);
  snew(stot,nprocs);
  init_nrnb(av);
  for(i=0; (i<nprocs); i++) {
    add_nrnb(av,av,&(nrnb[i]));
    for(j=0; (j<NFORCE_INDEX); j++) {
      ftot[i]+=nrnb[i].n[force_index[j]]*cost_nrnb(force_index[j]);
    }
    for(j=0; (j<NSHAKE_INDEX); j++) {
      stot[i]+=nrnb[i].n[shake_index[j]]*cost_nrnb(shake_index[j]);
    }
  }    
  for(j=0; (j<eNRNB); j++)
    av->n[j]=av->n[j]/(double)nprocs;
  
  fprintf(log,"\nDetailed load balancing info in percentage of average\n");
  
  fprintf(log,"Type        CPU:");
  for(i=0; (i<nprocs); i++)
    fprintf(log,"%3d ",i);
  fprintf(log,"Scaling\n");
  fprintf(log,"----------------");
  for(i=0; (i<nprocs); i++)
    fprintf(log,"----");
  fprintf(log,"-------\n");
  
  for(j=0; (j<eNRNB); j++) {
    unb=100.0;
    if (av->n[j] > 0) {
      fprintf(log,"%15s:",nrnb_str(j));
      for(i=0; (i<nprocs); i++) {
	dperc=(100.0*nrnb[i].n[j])/av->n[j];
	unb=max(unb,dperc);
	perc=dperc;
	fprintf(log,"%3d ",perc);
      }
      if (unb > 0) {
	perc=10000.0/unb;
	fprintf(log,"%6d%%\n",perc);
      }
      else
	fprintf(log,"\n");
    }   
  }
  fav=sav=0;
  for(i=0; (i<nprocs); i++) {
    fav+=ftot[i];
    sav+=stot[i];
  }
  uf=pr_av(log,nprocs,fav,ftot,"Total Force");
  us=pr_av(log,nprocs,sav,stot,"Total Shake");
  
  unb=(uf*fav+us*sav)/(fav+sav);
  if (unb > 0) {
    unb=10000.0/unb;
    fprintf(log,"\nTotal Scaling: %.0f%% of max performance\n\n",unb);
  }
}

