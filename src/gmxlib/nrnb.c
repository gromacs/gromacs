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
  { "LJ",                              31 }, /* inl0100 */
  { "LJ(S)",                           31 }, /* inl0110 */
  { "Buckingham",                      36 }, /* inl0200 */
  { "Buckingham(S)",                   36 }, /* inl0210 */
  { "LJ(T)",                           52 }, /* inl0300 */
  { "LJ(T)(S)",                        52 }, /* inl0310 */
  { "FreeEner LJ(T)",                  65 }, /* inl0301 */
  { "Softcore LJ(T)",                 118 }, /* inl0302 */
  { "Buckingham(T)",                   57 }, /* inl0400 */
  { "Buckingham(T)(S)",                57 }, /* inl0410 */
  { "FreeEner Bham(T)",                89 }, /* inl0401 */
  { "Softcore Bham(T)",               128 }, /* inl0402 */
  { "Coulomb",                         27 }, /* inl1000 */
  { "Coulomb(S)",                      37 }, /* inl1010 */
  { "Coulomb(W)",                      81 }, /* inl1020 */
  { "Coulomb(WW)",                    234 }, /* inl1030 */
  { "LJ + Coulomb",                    38 }, /* inl1100 */
  { "LJ + Coul(S)",                    41 }, /* inl1110 */
  { "LJ + Coul(W)",                    92 }, /* inl1120 */
  { "LJ + Coul(WW)",                  245 }, /* inl1130 */
  { "Buckingham + Coul",               40 }, /* inl1200 */
  { "Bham + Coul(S)",                  42 }, /* inl1210 */
  { "Bham + Coul(W)",                  94 }, /* inl1220 */
  { "Bham + Coul(WW)",                247 }, /* inl1230 */
  { "LJ(T) + Coul ",                   58 }, /* inl1300 */
  { "LJ(T) + Coul(S)",                 48 }, /* inl1310 */
  { "LJ(T) + Coul(W)",                112 }, /* inl1320 */
  { "LJ(T) + Coul(WW)",               265 }, /* inl1330 */
  { "Bham(T) + Coul ",                 63 }, /* inl1400 */
  { "Bham(T) + Coul(S)",               49 }, /* inl1410 */
  { "Bham(T) + Coul(W)",              117 }, /* inl1420 */
  { "Bham(T) + Coul(WW)",             270 }, /* inl1430 */
  { "RF Coul",                         32 }, /* inl2000 */
  { "RF Coul(S)",                      42 }, /* inl2010 */
  { "RF Coul(W)",                      96 }, /* inl2020 */
  { "RF Coul(WW)",                    279 }, /* inl2030 */
  { "LJ + RF Coul",                    43 }, /* inl2100 */
  { "LJ + RF Coul(S)",                 46 }, /* inl2110 */
  { "LJ + RF Coul(W)",                107 }, /* inl2120 */
  { "LJ + RF Coul(WW)",               290 }, /* inl2130 */
  { "Bham + RF Coul",                  45 }, /* inl2200 */
  { "Bham + RF Coul(S)",               47 }, /* inl2210 */
  { "Bham + RF Coul(W)",              109 }, /* inl2220 */
  { "Bham + RF Coul(WW)",             292 }, /* inl2230 */
  { "LJ(T) + RF Coul",                 63 }, /* inl2300 */
  { "LJ(T) + RF Coul(S)",              53 }, /* inl2310 */
  { "LJ(T) + RF Coul(W)",             127 }, /* inl2320 */
  { "LJ(T) + RF Coul(WW)",            310 }, /* inl2330 */
  { "Bham(T) + RF",                    68 }, /* inl2400 */
  { "Bham(T) + RF(S)",                 54 }, /* inl2410 */
  { "Bham(T) + RF(W)",                132 }, /* inl2420 */
  { "Bham(T) + RF(WW)",               310 }, /* inl2430 */
  { "Coulomb(T)",                      41 }, /* inl3000 */
  { "FreeEner Coul(T)",                54 }, /* inl3001 */
  { "Softcore Coul(T)",                99 }, /* inl3002 */
  { "Coulomb(T)(S)",                   51 }, /* inl3010 */
  { "Coulomb(T)(W)",                  123 }, /* inl3020 */
  { "Coulomb(T)(WW)",                 360 }, /* inl3030 */
  { "LJ + Coulomb(T)",                 54 }, /* inl3100 */
  { "LJ + Coulomb(T)(S)",              56 }, /* inl3110 */
  { "LJ + Coulomb(T)(W)",             136 }, /* inl3120 */
  { "LJ + Coulomb(T)(WW)",            373 }, /* inl3130 */
  { "Bham + Coul(T)",                  55 }, /* inl3200 */
  { "Bham + Coul(T)(S)",               56 }, /* inl3210 */
  { "Bham + Coul(T)(W)",              137 }, /* inl3220 */
  { "Bham + Coul(T)(WW)",             374 }, /* inl3230 */
  { "LJ(T) + Coul(T)",                 67 }, /* inl3300 */
  { "Free LJ(T)+Coul(T)",              92 }, /* inl3301 */
  { "SC LJ(T)+Coul(T)",               151 }, /* inl3302 */
  { "LJ(T) + Coul(T)(S)",              60 }, /* inl3310 */
  { "LJ(T) + Coul(T)(W)",             149 }, /* inl3320 */
  { "LJ(T) + Coul(T)(WW)",            386 }, /* inl3330 */
  { "LJ(T) + Coul(T)",                 71 }, /* inl3400 */
  { "Free Bham(T)+Coul(T)",           116 }, /* inl3401 */
  { "SC Bham(T)+Coul(T)",             161 }, /* inl3402 */
  { "Bham(T) + Coul(T)(S)",            61 }, /* inl3410 */
  { "Bham(T) + Coul(T,W)",            153 }, /* inl3420 */
  { "Bham(T) + Coul(T,WW)",           390 }, /* inl3430 */
  { "Innerloop-Iatom",                 10 },
  { "Calc Weights",                    36 },
  { "Spread Q",                         6 },
  { "Spread Q Bspline",                 2 }, 
  { "Gather F",                        23 },
  { "Gather F Bspline",                12 }, 
  { "3D-FFT",                           8 },
  { "Convolution",                      4 },
  { "Solve PME",                       64 },
  { "NS-Pairs",                        21 },
  { "Reset In Box",                     9 },
  { "Shift-X",                          6 },
  { "CG-CoM",                          29 },
  { "Sum Forces",                       1 },
  { "Bonds",                           43 },
  { "G96Bonds",                        40 },
  { "Angles",                         163 },
  { "G96Angles",                      150 },
  { "Propers",                        229 },
  { "Impropers",                      208 },
  { "RB-Dihedrals",                   247 },
  { "Dist. Restr.",                   200 },
  { "Pos. Restr.",                     50 },
  { "Angle Restr.",                   191 },
  { "Angle Restr. Z",                 164 },
  { "Morse Potent.",                   58 },
  { "Cubic Bonds",                     54 },
  { "Water Pol.",                      62 },
  { "Virial",                          18 },
  { "Update",                          31 },
  { "Ext.ens. Update",                 54 },
  { "Stop-CM",                         10 },
  { "P-Coupling",                       6 },
  { "Calc-Ekin",                       27 },
  { "Lincs",                           60 },
  { "Lincs-Mat",                        4 },
  { "Shake",                           30 },
  { "Shake-V",                         15 },
  { "Shake-Init",                      10 },
  { "Shake-Vir",                       18 },
  { "Settle",                         323 },
  { "PShake-InitLD",                   59 },    
  { "PShake-InitMD",                   65 },   
  { "PShake",                           7 },
  { "Dummy2",                          17 },
  { "Dummy3",                          28 },
  { "Dummy3fd",                        95 },
  { "Dummy3fad",                      176 },
  { "Dummy3out",                       87 },
  { "Dummy4fd",                       110 } 
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
  printf("nrnb %15s(%2d) incremented with %8d from file %s line %d\n",
	  nbdata[enr].name,enr,inc,file,line);
#endif
}

void print_perf(FILE *out,double nodetime,double realtime,real runtime,
		t_nrnb *nrnb,int nprocs)
{
  int    i;
  double nbfs,mni,frac,tfrac,mflop,tflop;
  
  if (nodetime == 0.0) {
    fprintf(out,"nodetime = 0! Infinite Giga flopses!\n");
  }
  
  nbfs=0.0;
  for(i=0; (i<eNR_INLOOP); i++) {
    if (strstr(nbdata[i].name,"(WW)") != NULL)
      nbfs += 9e-6*nrnb->n[i];
    else if (strstr(nbdata[i].name,"(W)") != NULL)
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
    nodetime = realtime;
    fprintf(out,"\tBased on real time for parallel computer.\n");
  }

  fprintf(out,"   RF=Reaction-field  Free=Free Energy  SC=Softcore\n");
  fprintf(out,"   T=Tabulated        S=Solvent         W=Water     WW=Water-Water\n\n");
  fprintf(out,"%21s  %12s  %12s  %8s\n",
	  "Computing:","M-Number","M-Flop's","% Flop's");
  mflop=0.0;
  tfrac=0.0;
  for(i=0; (i<eNRNB); i++) {
    mni    = 1e-6*nrnb->n[i];
    mflop += mni*nbdata[i].flop;
    frac   = 100.0*mni*nbdata[i].flop/tflop;
    tfrac += frac;
    if (mni != 0)
      fprintf(out,"%21s  %12.6f  %12.6f  %6.1f\n",
	      nbdata[i].name,mni,mni*nbdata[i].flop,frac);
  }
  fprintf(out,"%15s  %12s  %12.5f  %6.1f\n\n",
	  "Total","",mflop,tfrac);
  if (nodetime > 0 && realtime > 0) {
    fprintf(out,"%12s %10s %10s %8s\n","","NODE (s)","Real (s)","(%)");
    fprintf(out,"%12s %10.3f %10.3f %8.1f\n","Time:",
	    nodetime, realtime, 100.0*nodetime/realtime);
    if (nodetime > 60) {
      fprintf(out,"%12s %10s","","");
      pr_difftime(out,nodetime);
    }
    if (runtime>0) { /* runtime=0 means calc energies only */
      fprintf(out,"%12s %10s %10s %10s %10s\n",
	      "","(Mnbf/s)","(MFlops)","(ps/NODE hour)","(NODE hour/ns)");
      fprintf(out,"%12s %10.3f %10.3f %10.3f %10.3f\n","Performance:",
	      nbfs/nodetime,mflop/nodetime,
	      runtime*3600/nodetime,1000*nodetime/(3600*runtime));
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
  
  fprintf(log,"Type        NODE:");
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

