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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_nrnb_c = "$Id$";

#include <string.h>
#include "sysstuff.h"
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
} t_nrnb_data;

static t_nrnb_data nbdata[eNRNB] = {
  { "LJ+Coulomb",      38 },
  { "Coulomb",         27 },
  { "LJC-RF",          43 },
  { "Coul-RF",         32 },
  { "Buckingham",      44 },
  { "Buck.+RF",        49 },
  { "Table-Coul",      47 },
  { "Table-LJC",       73 },
  { "Table-BHAM",      78 },
  { "Table-BHAM-H2O", 171 },
  { "LJ+Coulomb-H2O",  91 },
  { "Coulomb-H2O",     80 },
  { "LJC-RF-H2O",     106 },
  { "Coul-RF-H2O",     95 },
  { "Buckingham-H2O",  97 },
  { "Buck.+RF-H2O",   112 },
  { "Table-Coul-H2O", 140 },
  { "Table-LJC-H2O",  166 },
  { "LJC-FreeEner",   101 },
  { "BHAM-FreeEner",  106 },
  { "LJC-Ewald",      150 }, /* not correct yet, but not used either */
  { "Coul-Ewald",     150 }, /* not correct yet, but not used either */
  { "BHAM-Ewald",     150 }, /* not correct yet, but not used either */
  { "LJC-Ewald-H2O",  150 }, /* not correct yet, but not used either */
  { "Coul-Ewald-H2O", 150 }, /* not correct yet, but not used either */
  { "BHAM-Ewald-H2O", 150 }, /* not correct yet, but not used either */
  { "Innerloop-Iatom", 10 },
  { "Calc Weights",    36 },
  { "Spread Q",         6 },
  { "Spread Q Bspline", 2 }, /* a first guess */
  { "Gather F",        23 },
  { "Gather F Bspline",12 }, /* a first guess */
  { "3D-FFT",           8 },
  { "Convolution",      4 },
  { "Solve PME",       64 }, /* a first guess */
  { "NS-Pairs",        21 },
  { "Reset In Box",     9 },
  { "Shift-X",          6 },
  { "CG-CoM",          29 },
  { "Sum Forces",       1 },
  { "Bonds",           43 },
  { "G96Bonds",        40 },
  { "Angles",         163 },
  { "G96Angles",      150 },
  { "Propers",        229 },
  { "Impropers",      208 },
  { "RB-Dihedrals",   247 },
  { "Dist. Restr.",   200 },
  { "Pos. Restr.",     50 },
  { "Angle Restr.",   191 },
  { "Angle Restr. Z", 164 },
  { "Morse Potent.",   58 },
  { "Cubic Bonds",     54 },
  { "Water Pol.",      62 },
  { "Virial",          18 },
  { "Update",          31 },
  { "Stop-CM",         10 },
  { "P-Coupling",       3 },
  { "Calc-Ekin",       27 },
  { "Lincs",           60 },
  { "Lincs-Mat",        4 },
  { "Shake",           30 },
  { "Shake-V",         15 },
  { "Shake-Init",      10 },
  { "Shake-Vir",       18 },
  { "Settle",         323 },
  { "PShake-InitLD",   59 },    
  { "PShake-InitMD",   65 },   
  { "PShake",           7 },
  { "Dummy2",          17 },
  { "Dummy3",          28 },
  { "Dummy3fd",        95 },
  { "Dummy3fad",      176 },
  { "Dummy3out",       87 },
  { "Dummy4fd",       110 }
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

void print_perf(FILE *out,double cputime,double realtime,real runtime,
		t_nrnb *nrnb,int nprocs)
{
  int    i;
  double nbfs,mni,frac,tfrac,mflop,tflop;
  
  if (cputime == 0.0) {
    fprintf(out,"cputime = 0! Infinite Giga flopses!\n");
  }
  
  nbfs=0.0;
  for(i=0; (i<eNR_INLOOP); i++) {
    if (strstr(nbdata[i].name,"H2O") != NULL)
      nbfs += 3e-6*nrnb->n[i];
    else
      nbfs += 1e-6*nrnb->n[i];
  }
  tflop=0;
  for(i=0; (i<eNRNB); i++) 
    tflop+=1e-6*nrnb->n[i]*nbdata[i].flop;
  
  if (tflop == 0) {
    fprintf(out,"No MEGA Flopsen this time\n");
    return;
  }
  fprintf(out,"\n\tM E G A - F L O P S   A C C O U N T I N G\n\n");
  if (nprocs > 1) {
    cputime = realtime;
    fprintf(out,"\tBased on real time for parallel computer.\n");
  }
  fprintf(out,"%15s  %12s  %12s  %8s\n",
	  "Computing:","M-Number","M-Flop's","% Flop's");
  mflop=0.0;
  tfrac=0.0;
  for(i=0; (i<eNRNB); i++) {
    mni    = 1e-6*nrnb->n[i];
    mflop += mni*nbdata[i].flop;
    frac   = 100.0*mni*nbdata[i].flop/tflop;
    tfrac += frac;
    if (mni != 0)
      fprintf(out,"%15s  %12.6f  %12.6f  %6.1f\n",
	      nbdata[i].name,mni,mni*nbdata[i].flop,frac);
  }
  fprintf(out,"%15s  %12s  %12.5f  %6.1f\n\n",
	  "Total","",mflop,tfrac);
  if (cputime > 0 && realtime > 0) {
    fprintf(out,"%12s %10s %10s %8s\n","","CPU (s)","Real (s)","(%)");
    fprintf(out,"%12s %10.3f %10.3f %8.1f\n","Time:",
	    cputime, realtime, 100.0*cputime/realtime);
    if (cputime > 60) {
      fprintf(out,"%12s %10s","","");
      pr_difftime(out,cputime);
    }
    if (runtime>0) { /* runtime=0 means calc energies only */
      fprintf(out,"%12s %10s %10s %10s %10s\n",
	      "","(Mnbf/s)","(MFlops)","(ps/CPU hour)","(CPU hour/ns)");
      fprintf(out,"%12s %10.3f %10.3f %10.3f %10.3f\n","Performance:",
	      nbfs/cputime,mflop/cputime,
	      runtime*3600/cputime,1000*cputime/(3600*runtime));
    }
  }
}

int cost_nrnb(int enr)
{
  return nbdata[enr].flop;
}

char *nrnb_str(int enr)
{
  return nbdata[enr].name;
}

static int    force_index[]={ 
  eNR_BONDS,  eNR_ANGLES,  eNR_PROPER, eNR_IMPROPER, 
  eNR_RB,     eNR_DISRES,  eNR_POSRES,
  eNR_NS,     eNR_INL_IATOM
};
#define NFORCE_INDEX asize(force_index)

static int    shake_index[]={ 
  eNR_SHAKE,     eNR_SHAKE_RIJ, eNR_SETTLE,       eNR_UPDATE,       eNR_PCOUPL,
  eNR_SHAKE_VIR, eNR_SHAKE_V,   eNR_PSHAKEINITLD, eNR_PSHAKEINITMD, eNR_PSHAKE
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
    /* Cost due to forces */
    for(j=0; (j<eNR_INLOOP); j++)
      ftot[i]+=nrnb[i].n[j]*cost_nrnb(j);
    for(j=0; (j<NFORCE_INDEX); j++) 
      ftot[i]+=nrnb[i].n[force_index[j]]*cost_nrnb(force_index[j]);
    /* Due to shake */
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

